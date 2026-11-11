"""
GPU-accelerated, parallel metabolic control analysis over many parameter sets.

    import tellurium as te, numpy as np
    te.gpu(True)                       # CUDA (NVIDIA), MPS (Apple Silicon), or CPU fallback
    r = te.loada('...')
    res = te.controlCoefficients(r, {'k1': np.linspace(0.1, 2, 10000)})
    res['fcc'].shape                   # (10000, nReactions, nReactions)

How the work is split:
  * steady states + unscaled elasticities: the Rust engine ``roadrunner_rs`` (rrs/) on all
    cores when it is built and supports the model, else roadrunner in parallel worker processes
  * C = -L (N_R E L)^-1 N_R,  F = I + E C,  and scaling: one batched solve on the GPU

Scripts (not notebooks) that use more than one process need the usual
``if __name__ == '__main__':`` guard, since workers are spawned.
"""
from __future__ import absolute_import
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool

import numpy as np

_device = None  # torch.device when GPU is on, None means numpy on CPU


def gpu(enable=True):
    """Turn GPU acceleration for :func:`controlCoefficients` on or off.

    :param enable: True/False, "true"/"false", or a device name ("cuda", "mps", "cpu").
    :returns: name of the device now in use ("cuda", "mps" or "cpu").

    Needs PyTorch (``pip install torch``). Falls back to CPU with a warning when
    no GPU is found. Apple MPS only supports float32, so results there agree with
    the CPU to ~1e-5 relative, not 1e-12.
    """
    global _device
    want = str(enable).strip().lower()
    if want in ('false', '0', 'no', 'off', 'none', 'cpu'):
        _device = None
        return 'cpu'
    try:
        import torch
    except ImportError:
        warnings.warn("te.gpu: PyTorch is not installed (pip install torch); using CPU.")
        _device = None
        return 'cpu'
    found = 'cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else None
    if want in ('cuda', 'mps') and want != found:
        found = None
    if found is None:
        warnings.warn("te.gpu: no usable GPU found; using CPU.")
        _device = None
        return 'cpu'
    _device = torch.device(found)
    return found


def getGPU():
    """Name of the device used by :func:`controlCoefficients` ("cuda", "mps" or "cpu")."""
    return 'cpu' if _device is None else _device.type


# ---- CPU worker side: steady state + elasticities per parameter set ----

_rr = None


def _init_worker(sbml):
    global _rr
    from roadrunner.roadrunner import RoadRunner  # plain libroadrunner, not tellurium's Rust-backed subclass
    _rr = RoadRunner(sbml)
    # always on: with it off, every steadyState() rebuilds + re-JITs the model twice (2.4x slower)
    _rr.conservedMoietyAnalysis = True


def _sample(params):
    _rr.resetAll()
    for k, v in params.items():
        _rr[k] = float(v)  # numpy ints are rejected by roadrunner
    try:
        _rr.steadyState()
        return (np.array(_rr.getUnscaledElasticityMatrix()),
                np.array(_rr.getReactionRates()),
                np.array(_rr.model.getFloatingSpeciesConcentrations()))
    except Exception:
        return None  # no steady state for this parameter set -> NaN row


# ---- batched linear algebra (GPU or numpy) ----

def _mca(NR, L, E, v, S):
    if _device is not None:
        import torch
        dt = torch.float32 if _device.type == 'mps' else torch.float64
        NR, L, E, v, S = (torch.as_tensor(a, dtype=dt, device=_device) for a in (NR, L, E, v, S))
        solve, eye = torch.linalg.solve, torch.eye(E.shape[1], dtype=dt, device=_device)
    else:
        solve, eye = np.linalg.solve, np.eye(E.shape[1])
    C = -L @ solve(NR @ E @ L, NR)                     # (B, species, reactions)
    F = eye + E @ C                                    # (B, reactions, reactions)
    out = (C * v[:, None, :] / S[:, :, None],          # scaled CCC
           F * v[:, None, :] / v[:, :, None],          # scaled FCC
           E * S[:, None, :] / v[:, :, None],          # scaled elasticities
           C, F)
    if _device is not None:
        out = tuple(a.cpu().double().numpy() for a in out)
    return out


def controlCoefficients(r, samples=None, processes=None):
    """Concentration/flux control coefficients and elasticities for many parameter sets at once.

    :param r: roadrunner model (e.g. from ``te.loada``). It is not modified.
    :param samples: parameter sets to evaluate. Either a dict of equal-length
        sequences ``{'k1': [...], 'k2': [...]}`` (zipped, not a grid), a list of
        dicts, or a pandas DataFrame. None evaluates the model's current values once.
    :param processes: CPU threads (Rust engine) or worker processes (roadrunner fallback)
        for the steady states. Default: all cores (Rust); roadrunner: all cores at 64+ samples, else 1.
    :returns: dict of numpy arrays, first axis = sample:
        ``ccc``/``fcc``/``elasticities`` (scaled), ``uccc``/``ufcc`` (unscaled),
        plus ``species``, ``reactions`` id lists and ``ok`` (False where no
        steady state was found; those rows are NaN).

    Use ``te.gpu(True)`` first to run the matrix algebra on CUDA or Apple MPS.
    """
    if samples is None:
        samples = [{}]
    elif hasattr(samples, 'to_dict'):
        samples = samples.to_dict('records')
    elif isinstance(samples, dict):
        keys = list(samples)
        samples = [dict(zip(keys, row)) for row in zip(*(np.atleast_1d(samples[k]) for k in keys))]
    samples = list(samples)
    n = len(samples)

    keys = list(samples[0]) if n else []
    m = None
    if all(list(p) == keys for p in samples):
        try:
            from .rrs_bridge import from_roadrunner
            m = from_roadrunner(r)
        except ImportError:
            pass  # roadrunner_rs not built: roadrunner path
        except NotImplementedError as e:
            warnings.warn("te.controlCoefficients: Rust engine can't run this model (%s); using roadrunner." % e)
    if m is None:
        return _cc_roadrunner(r, samples, processes)

    # Rust engine: steady states + exact elasticities on all cores, no worker processes
    vals = np.array([[float(p[k]) for k in keys] for p in samples]).reshape(n, len(keys))
    # ponytail: GPU ignored here; the per-sample algebra is ~5 us on CPU, a GPU copy costs more
    try:
        res = m.mca_batch(keys, vals, threads=processes)
    except ValueError:  # rate rules (or unknown key: roadrunner then raises its own error)
        return _cc_roadrunner(r, samples, processes)
    out = [res[k] for k in _KEYS]
    ok = res['ok']
    bad = np.flatnonzero(~ok)
    if bad.size:  # retry Rust failures with roadrunner's solvers
        rr = _cc_roadrunner(r, [samples[i] for i in bad], processes)
        for dst, k in zip(out, _KEYS):
            dst[bad] = rr[k]
        ok = ok.copy()
        ok[bad] = rr['ok']
    return dict(zip(_KEYS, out), species=list(m.species), reactions=list(m.reactions), ok=ok)


_KEYS = ('ccc', 'fcc', 'elasticities', 'uccc', 'ufcc')


def _cc_roadrunner(r, samples, processes):
    n = len(samples)
    sbml = r.getCurrentSBML()
    # ponytail: process pool only pays off once per-worker model compile (~0.1-0.5 s) is amortized
    processes = processes or (os.cpu_count() if n >= 64 else 1)
    _init_worker(sbml)  # local copy: serial runs + structural matrices below
    res = None
    if processes > 1:
        try:
            with ProcessPoolExecutor(processes, initializer=_init_worker, initargs=(sbml,)) as ex:
                res = list(ex.map(_sample, samples, chunksize=max(1, n // (4 * processes))))
        except BrokenProcessPool:
            warnings.warn("te.controlCoefficients: worker processes failed to start (a script needs an "
                          "`if __name__ == '__main__':` guard); falling back to 1 process.")
    if res is None:
        res = [_sample(p) for p in samples]

    # Structure is parameter independent: take it once. Align link-matrix rows to elasticity columns.
    # Moiety analysis reorders species; report them in the model's own order.
    Em = _rr.getUnscaledElasticityMatrix()
    species, reactions = list(r.model.getFloatingSpeciesIds()), list(Em.rownames)
    col = [list(Em.colnames).index(s) for s in species]
    Lm = _rr.getLinkMatrix()
    L = np.array(Lm)[[list(Lm.rownames).index(s) for s in species]]
    NR = np.array(_rr.getReducedStoichiometryMatrix())

    ok = np.array([x is not None for x in res])
    nS, nR = len(species), len(reactions)
    shapes = [(nS, nR), (nR, nR), (nR, nS), (nS, nR), (nR, nR)]
    out = [np.full((n,) + s, np.nan) for s in shapes]
    if ok.any():
        good = [x for x in res if x is not None]
        E, v, S = (np.stack(a) for a in zip(*good))
        E, S = E[:, :, col], S[:, col]
        for dst, src in zip(out, _mca(NR, L, E, v, S)):
            dst[ok] = src

    return dict(zip(_KEYS, out), species=species, reactions=reactions, ok=ok)
