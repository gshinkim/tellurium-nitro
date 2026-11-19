"""
Extending the RoadRunner object with additional fields.
"""

from __future__ import print_function, division, absolute_import
import os
import roadrunner
import warnings
import copy
import inspect
import numpy as np

# ---------------------------------------------------------------------
# Turbo: route roadrunner's expensive calls to the Rust engine roadrunner_rs
# ---------------------------------------------------------------------
try:
    import roadrunner_rs  # noqa: F401
    _has_rs = True
    _turbo = os.environ.get('TELLURIUM_ENGINE', '').strip().lower() != 'roadrunner'
except ImportError:
    _has_rs = _turbo = False
# strict: never hand a model or call to libroadrunner; raise NotSupportedByRust instead (te.turbo('strict'))
_strict = _has_rs and os.environ.get('TELLURIUM_ENGINE', '').strip().lower() == 'strict'

_NamedArray = roadrunner._roadrunner.NamedArray


class _Opts:  # the SimulateOptions fields _rs_simulate reads, for the simulate fast path
    times, output_file, reset_model = (), '', False

    def __init__(self, start, end, steps):
        self.start, self.end, self.steps = start, end, steps
# integrator settings the Rust engine honours only at these values (others fall back to roadrunner).
# cvode: tolerances are passed through; stiff/BDF/Adams orders/maximum_num_steps don't change the
# solution within tolerance (Rodas4 is stiff-stable and has its own step cap).
_RS_INTEGRATORS = {
    'cvode': {'maximum_time_step': 0.0, 'minimum_time_step': 0.0, 'initial_time_step': 0.0,
              'multiple_steps': False, 'variable_step_size': False},
    'gillespie': {'minimum_time_step': 0.0, 'maximum_time_step': 0.0, 'nonnegative': False, 'maximum_num_steps': 0},
}
_RS_INTEGRATORS = {k: (tuple(v), tuple(v.values())) for k, v in _RS_INTEGRATORS.items()}  # (keys, values)


def turbo(enable=None):
    """Use the Rust engine (roadrunner_rs) for simulate/steadyState/MCA when it supports the model.

    :param enable: True/False to switch globally, 'strict' for Rust only (anything the Rust engine
        can't do raises NotSupportedByRust instead of going to roadrunner); None only queries.
    :returns: whether turbo is on. On by default when roadrunner_rs is installed;
        env ``TELLURIUM_ENGINE=roadrunner`` starts it off, ``TELLURIUM_ENGINE=strict`` starts it strict.
        Otherwise unsupported models, integrators or settings fall back to roadrunner.
    """
    global _turbo, _strict
    if enable is not None:
        _strict = enable == 'strict' and _has_rs
        if enable:
            try:
                import roadrunner_rs  # noqa: F401
            except ImportError:
                warnings.warn("te.turbo: roadrunner_rs is not installed; using roadrunner.")
                enable = False
        _turbo = bool(enable)
    return _turbo


def _named(a, rows=None, cols=None):
    a = np.ascontiguousarray(a, dtype=float).view(_NamedArray)
    if rows is not None:
        a.rownames = list(rows)
    if cols is not None:
        a.colnames = list(cols)
    return a


class _Fallback(Exception):
    """Raised inside a Rust-engine path: the call goes to roadrunner instead; the message is the reason."""


class NotSupportedByRust(RuntimeError):
    """te.turbo('strict'): the Rust engine can't do this and roadrunner fallback is off."""


def _no_fallback(reason):
    if _strict:
        raise NotSupportedByRust(str(reason))


class _Cols(object):
    """Plan turning the Rust model's output columns into roadrunner's selections (one per selection list)."""

    def __init__(self, sel, outs, ores, I, mul, div, lin, init):
        self.sel, self.outs, self.ores = sel, outs, ores
        self.I, self.mul, self.div, self.lin, self.init = I, mul, div, lin, init
        self.ident = I == list(range(len(outs))) and not (mul or div or lin or init)

    def init_values(self, r):
        """init(...) selections: roadrunner's values now (taken before a run: they need roadrunner's model)."""
        return [roadrunner.RoadRunner.getValue(r, s) for _, s in self.init]

    def apply(self, Y, pre=()):
        """Rust output rows Y (columns = self.outs) -> rows in selection order."""
        if self.ident:
            return Y
        out = Y[:, self.I]
        for ops, f in ((self.mul, np.multiply), (self.div, np.divide)):
            if ops:
                k, j = zip(*ops)
                out[:, list(k)] = f(out[:, list(k)], Y[:, list(j)])
        for k, j, c in self.lin:
            out[:, k] = Y[:, j] @ c
        for (k, _), x in zip(self.init, pre):
            out[:, k] = x
        return out


class _RsModel(object):
    """The Rust model ``m`` of one roadrunner model (built by rrs_bridge from ``info``), the index maps
    between roadrunner's ids and it, and the running simulation (event state, solver counters).
    ``ids`` = (floating, boundary, parameter, compartment, reaction ids); ``vols``: compartment
    volumes; ``nss``: no steady state (time-dependent rates or events)."""
    rd = L = ind = tk = ick = None  # filled on first use

    def __init__(self, m, info, ids, vols, nss, extra=False):
        ids, bids, gids, cids, rids = ids
        states, outs = list(m.state_ids), set(m.output_ids)
        if not set(ids) <= outs or sorted(rids) != sorted(m.reactions):
            raise _Fallback("roadrunner's species/reactions are not all in the Rust model")

        def known(k):  # Rust parameter or state (not a rule output / moiety total)
            try:
                m.resolve_keys([k])
                return True
            except ValueError:
                return False
        pidx, cidx, bidx, sidx = ([i for i, k in enumerate(l) if known(k)] for l in (gids, cids, bids, ids))
        keys = [gids[i] for i in pidx] + [cids[i] for i in cidx] + [bids[i] for i in bidx]
        # extra (lazy models only): Rust states in none of roadrunner's groups, e.g. a species reference
        # with a rate rule; they sit before the floating species and have no roadrunner array (xidx)
        self.n_extra = len(keys)
        keys += [k for k in states if extra and k not in keys and k not in ids]
        self.n_extra = len(keys) - self.n_extra
        keys += [ids[i] for i in sidx]
        comp_of, general = info['comp_of'], info['general']
        rx = [list(m.reactions).index(j) for j in rids]
        NR, Lc, _ = m.conservation()
        Ngen = np.asarray(Lc) @ np.asarray(NR)  # amount stoichiometry, states x (reactions + rate rules)
        self.m, self.session = m, roadrunner_rs.Session(m)
        self.kres, self.ids, self.rids, self.gids, self.cids, self.bids = m.resolve_keys(keys), ids, rids, gids, cids, bids
        self.pidx, self.cidx, self.bidx, self.sidx, self.keys = pidx, cidx, bidx, sidx, keys
        self.kpos, self.off = {k: i for i, k in enumerate(keys)}, len(keys) - len(sidx)
        self.bamt = np.array([bids[i] in info['amount_params'] for i in bidx], dtype=bool)
        self.amount_params, self.comp_of, self.general, self.outs, self.states = info['amount_params'], comp_of, general, outs, states
        self.rx, self.cols_cache, self.nss = rx, {}, nss
        self.sp = [states.index(s) if s in states else None for s in ids]
        self.asc = np.ones(len(states))
        self.V = np.array([vols[cids.index(comp_of[s])] for s in ids])
        # species volumes the Rust model bakes in (rebuild when they change)
        self.cconst = sorted({cids.index(comp_of[s]) for s in m.species if s not in general})
        self.spi = [i for i, j in enumerate(self.sp) if j is not None]
        self.spj = [self.sp[i] for i in self.spi]
        self.stats = None
        oid = {n: i for i, n in enumerate(m.output_ids)}
        nb = len(bids) if bidx else 0  # x() gathers [globals | volumes | boundary conc | boundary amounts | floating]
        o1, o2 = len(gids), len(gids) + len(cids)
        self.xidx = np.array(pidx + [o1 + i for i in cidx] + [o2 + i + nb * a for i, a in zip(bidx, self.bamt)] +
                             [o2 + 2 * nb + i for i in sidx], dtype=int)

        def plan(chg):  # what sync writes back: species, and parameters/compartments/boundary species in chg
            return [(f, idx, [names[i] for i in idx], np.array([oid[names[i]] for i in idx], dtype=int)) for f, idx, names in (
                ('setCompartmentVolumes', [i for i in cidx if cids[i] in chg], cids),  # volumes first: species keep conc
                ('setGlobalParameterValues', [i for i in pidx if gids[i] in chg], gids),
                ('amount', [i for i, a in zip(bidx, self.bamt) if bids[i] in chg and a], bids),
                ('conc', [i for i, a in zip(bidx, self.bamt) if bids[i] in chg and not a], bids),  # via setValue
                ('setFloatingSpeciesConcentrations', sidx, ids)) if idx]
        self.sync, self.sync_all = plan(info['changing']), plan(set(keys))
        # reaction-network views (MCA, Jacobian, elasticities) over the floating species; rate-rule
        # states are not floating species in roadrunner and are held fixed, as roadrunner does
        self.plain = all(s in m.species and s not in general for s in ids)
        self.mca = self.plain and not nss
        if self.plain:
            self.N = Ngen[self.sp][:, rx]

    def reset_events(self):
        """The next run starts with fresh event triggers (as after roadrunner's reset())."""
        self.session.reset_events()

    def times(self, start, end, steps):
        tk = (start, end, steps)
        if self.tk != tk:
            self.tk, self._times = tk, np.linspace(start, end, steps + 1)
        return self._times

    def atol_scale(self, ic):
        """Per-state atol multipliers from the floating species' initial concentrations ``ic``:
        roadrunner's CVODE multiplies a scalar atol by each species' initial amount (by its volume
        if that is zero), i.e. atol * |initial concentration| (atol if zero). Rate-rule states keep atol."""
        k = ic.tobytes()
        if self.ick != k:
            a = np.abs(ic[self.spi])
            self.ick, self.asc[self.spj] = k, np.where(a < 1e-30, 1.0, a)
        return self.asc

    def advance(self, x, times, ores, rtol, atol, asc):
        """Integrate from state x over times -> (outputs at times, final values of all outputs).
        Keeps the event state for the next call; self.stats = the solver's counters."""
        try:
            y, final = self.session.advance(self.kres, x, times, outputs=ores, rtol=rtol, atol=atol, atol_scale=asc)
        except RuntimeError:
            self.stats = self.session.stats
            raise _Fallback("Rust integration failed (atol %g and 1e-12)" % atol)
        self.stats = self.session.stats
        return np.asarray(y), final

    def cols(self, sel, lazy=False):
        """_Cols for the selection list, or None if Rust can't produce it (``lazy``: no roadrunner model,
        so no init(...) selections)."""
        key = tuple(sel)
        c = self.cols_cache.get(key, False)
        if c is False:
            c = self.cols_cache[key] = self._cols(list(sel))
        return None if c is None or (lazy and c.init) else c

    def _cols(self, sel):
        outs = []

        def o(n):
            if n not in outs:
                outs.append(n)
            return outs.index(n)
        sp = set(self.ids) | set(self.bids)
        I, mul, div, lin, init = [], [], [], [], []
        for k, s in enumerate(sel):
            if not isinstance(s, str):
                return None
            br = s[:1] == '[' and s[-1:] == ']'
            v = s[1:-1] if br else s
            if v in sp:
                if v not in self.outs:
                    return None
                I.append(o(v))
                amt = v in self.amount_params  # boundary species whose Rust value is an amount
                if br and amt:
                    div.append((k, o(self.comp_of[v])))
                elif not br and not amt:
                    mul.append((k, o(self.comp_of[v])))
            elif s in self.outs:  # time, parameters, compartments, reaction rates, rule variables
                I.append(o(s))
            elif s.startswith('init(') and s.endswith(')'):
                I.append(0)
                init.append((k, s))
            elif s.endswith("'") and self.plain and s[:-1].strip('[]') in self.ids:
                i = self.ids.index(s[:-1].strip('[]'))
                I.append(0)  # amount (or conc: constant compartment, plain) rate from the reaction rates
                lin.append((k, [o(j) for j in self.rids], self.N[i] / (self.V[i] if s[0] == '[' else 1.0)))
            else:
                return None  # ponytail: eigen(), cc(), ee() ... selections -> roadrunner
        return _Cols(sel, outs, self.m.resolve_outputs(outs), I, mul, div, lin, init)


# ---------------------------------------------------------------------
# Extended RoadRunner class
# ---------------------------------------------------------------------
class ExtendedRoadRunner(roadrunner.RoadRunner):

    def __new__(cls, *args, **kwargs):
        # RoadRunner(sbml) with turbo: a lazily compiled RoadRunner (see _lazy); otherwise a plain one
        r = cls._lazy(args[0]) if len(args) == 1 and not kwargs else None
        return r if r is not None else super(ExtendedRoadRunner, cls).__new__(cls)

    def __init__(self, *args, **kwargs):
        if '_lz' not in self.__dict__:  # lazy instances are already built by __new__
            super(ExtendedRoadRunner, self).__init__(*args, **kwargs)

    # ---------------------------------------------------------------------
    # Lazy roadrunner: with turbo, te.loada & co. build only the Rust model (~0.2 s at 1000 reactions)
    # and keep the model state here (lz['x'], lz['t']). roadrunner compiles the model (~10 s at 1000
    # reactions) on the first call the Rust engine can't answer (_load_rr), from the SBML + current state.
    # ---------------------------------------------------------------------
    @classmethod
    def _lazy(cls, sbml):
        """A lazily compiled RoadRunner when turbo is on and the Rust engine takes the model, else None."""
        if _turbo and isinstance(sbml, str) and not roadrunner.Config.getValue(
                roadrunner.Config.LOADSBMLOPTIONS_CONSERVED_MOIETIES):
            text, base = sbml, None
            if not text.lstrip().startswith('<') and os.path.isfile(text):
                base = os.path.dirname(os.path.abspath(text))
                with open(text, encoding='utf-8') as f:
                    text = f.read()
            r = object.__new__(cls)
            roadrunner.roadrunner.RoadRunner._swig_init(r)  # empty roadrunner (no _makeProperties: no model)
            try:
                if text.lstrip().startswith('<'):
                    r._lz_init(text, base)
                    return r
            except (NotImplementedError, ValueError, _Fallback) as e:
                _no_fallback('model not supported by the Rust engine: %s' % e)
        elif _turbo:
            _no_fallback('model needs roadrunner (conserved moiety analysis or not SBML text)')
        return None

    def _lz_init(self, sbml, base=None):
        from ..analysis.rrs_bridge import from_sbml, _names
        info = {}
        m = from_sbml(sbml, info, base)
        fids, bids, gids, cids, rids = ids = info['ids']
        doc = info['doc'].getModel()
        st = _RsModel(m, info, ids, [info['values'][c] for c in cids], 'symbols/time' in sbml or doc.getNumEvents() > 0,
                      extra=True)
        if not set(st.states) <= set(st.keys) or not set(sum(ids, [])) <= st.outs:
            raise NotImplementedError("not all of the model is in the Rust model")
        oid = {n: i for i, n in enumerate(m.output_ids)}
        x0 = np.asarray(m.outputs_at([], [], 0.0, outputs=st.keys), dtype=float)
        rr = info['rate_rule_ids']
        # reset(): floating species and rate-rule targets back to their initial values, an initial
        # assignment re-evaluated with the current parameters and initial species (as roadrunner)
        tgt = [k for k in st.keys if k in fids or k in rr]
        ok = not (set(cids) & (info['changing'] | set(info['ia'])))  # constant compartments
        plain = set(gids) - {x.getVariable() for x in doc.getListOfRules()}  # parameters without rules
        for k in tgt:
            if k in info['ia']:
                ok = ok and _names(info['ia'][k], set()) <= set(fids) | set(cids) | set(info['fdefs']) | plain
        sp = {s: info['doc'].getModel().getSpecies(s) for s in fids + bids}
        attr = {}
        for s in fids + bids:
            attr.update({s: '[%s]' % s, s + '_conc': '[%s]' % s, s + '_amt': s})
        attr.update({k: k for k in gids + cids + rids})
        srefs = [sr for rx in doc.getListOfReactions() for l in (rx.getListOfReactants(), rx.getListOfProducts()) for sr in l]
        self.__dict__['_lz'] = dict(
            sbml=sbml, info=info, st=st, x=x0.copy(), x0=x0, t=0.0, oid=oid, xo=np.array([oid[k] for k in st.keys], dtype=int),
            o=lambda l: np.array([oid[k] for k in l], dtype=int), attr=attr, reset_ok=ok,
            srefs=any(sr.isSetStoichiometryMath() for sr in srefs),  # these SBML views: roadrunner
            # named stoichiometries: (row, reaction column, sign, id); their current values update N
            nref=[(info['rows'].index(sr.getSpecies()), rids.index(sr.getParentSBMLObject().getParentSBMLObject().getId()),
                   -1 if sr.getParentSBMLObject().getElementName() == 'listOfReactants' else 1, sr.getId())
                  for sr in srefs if sr.isSetId() and sr.getSpecies() in info['rows']],
            tgt=np.array([st.kpos[k] for k in tgt], dtype=int), tgt_ids=tgt, model=_LazyModel(self),
            hosu={s for s in sp if sp[s].getHasOnlySubstanceUnits()},
            # roadrunner bakes these into the model: changing them needs roadrunner
            fixed=cids + [c for c in (doc.getConversionFactor(), *(sp[s].getConversionFactor() for s in sp)) if c],
            sel=['time'] + ['[%s]' % s for s in fids] + [v for v in rr if v not in fids], ssel=['[%s]' % s for s in fids])

    def _load_rr(self):
        """Compile the model in roadrunner (once) and carry the lazy state over; returns self."""
        lz = self.__dict__.get('_lz')
        if lz is None:
            return self
        _no_fallback('this call needs roadrunner')
        del self.__dict__['_lz']
        st = lz['st']
        roadrunner.roadrunner.RoadRunner.load(self, lz['sbml'])
        mdl = self.model
        now = self._rs_x(st, mdl)
        self._rs_sync(st, {k: v for k, v, w in zip(st.keys, lz['x'], now) if v != w}, mdl, plan='sync_all')
        mdl.setTime(lz['t'])
        roadrunner.roadrunner.RoadRunner._setSelections(self, lz['sel'])
        roadrunner.roadrunner.RoadRunner._setSteadyStateSelections(self, lz['ssel'])
        ids = [list(f()) for f in (mdl.getFloatingSpeciesIds, mdl.getBoundarySpeciesIds, mdl.getGlobalParameterIds,
                                   mdl.getCompartmentIds, mdl.getReactionIds)]
        if not st.n_extra and ids == [getattr(st, k) for k in ('ids', 'bids', 'gids', 'cids', 'rids')]:  # keep the Rust model
            cc = st.cconst
            self.__dict__['_rs_cache'] = ((self.conservedMoietyAnalysis, mdl.getCompartmentVolumes()[cc].tobytes()), st, cc)
        return self

    def _super(self):
        """roadrunner's own methods (compiling the model first if it is lazy)."""
        return super(ExtendedRoadRunner, self._load_rr())

    def _lz_out(self, lz, outputs=None):
        """Rust outputs (all, or the resolved `outputs`) at the current state and time."""
        st = lz['st']
        return np.asarray(st.m.outputs_at(st.kres, lz['x'], lz['t'], outputs=outputs))

    def _lz_get(self, lz, sel):
        cols = lz['st'].cols([sel], lazy=True) if isinstance(sel, str) else None
        if cols is None:
            return None
        return float(cols.apply(self._lz_out(lz, cols.ores)[None])[0, 0])

    def _lz_set(self, lz, sel, v):
        st = lz['st']
        br = isinstance(sel, str) and sel[:1] == '[' and sel[-1:] == ']'
        k = sel[1:-1] if br else sel
        if k == 'time' and not br:
            lz['t'] = float(v)
            return True
        j = st.kpos.get(k)
        sp = k in st.ids or k in st.bids
        if j is None or k in lz['fixed'] or (br and not sp):
            return False
        if sp:
            V = self._lz_out(lz)[lz['oid'][st.comp_of[k]]]
            c = float(v) if br else float(v) / V
            lz['x'][j] = c * V if k in st.amount_params else c
        else:
            lz['x'][j] = float(v)
        return True

    def _lz_init_values(self, lz):
        """Values reset() gives the floating species and rate-rule targets (lz['tgt_ids'])."""
        info, st, x0 = lz['info'], lz['st'], lz['x0']
        vals = dict(zip(lz['tgt_ids'], x0[lz['tgt']]))
        ia = {k: info['ia'][k] for k in lz['tgt_ids'] if k in info['ia']}
        if ia:
            from ..analysis.rrs_bridge import evaluate
            cur = dict(zip(st.keys, lz['x']))
            V = {c: cur[c] for c in st.cids if c in cur}

            def get(n):  # SBML units: species amounts when hasOnlySubstanceUnits
                if n in ia:
                    if n not in done:
                        done[n] = float(evaluate(ia[n], get, info['fdefs']))
                    v = done[n]
                    return v
                if n in st.ids:
                    return vals[n] * V[st.comp_of[n]] if n in lz['hosu'] else vals[n]
                return cur[n]
            done = {}
            for k in ia:
                v = get(k)
                vals[k] = v / V[st.comp_of[k]] if k in lz['hosu'] and k in st.ids else v
        return np.array([vals[k] for k in lz['tgt_ids']], dtype=float)

    def _lz_reset(self, lz, everything=False):
        lz['st'].reset_events()
        lz['t'] = 0.0
        if everything:
            lz['x'][:] = lz['x0']
        else:
            lz['x'][lz['tgt']] = self._lz_init_values(lz)

    def getModel(self):
        lz = self.__dict__.get('_lz')
        return roadrunner.roadrunner.RoadRunner.getModel(self) if lz is None else lz['model']
    getModel.__doc__ = roadrunner.roadrunner.RoadRunner.getModel.__doc__

    model = property(getModel)

    def _selections(name, base):
        def get(self):
            lz = self.__dict__.get('_lz')
            return list(lz[name]) if lz is not None else base.fget(self)

        def put(self, sel):
            lz = self.__dict__.get('_lz')
            if lz is not None:  # selections Rust computes, in roadrunner's spelling
                sel = list(sel)
                if lz['st'].cols(sel, lazy=True) is not None:
                    lz[name] = sel
                    return
                self._load_rr()
            base.fset(self, sel)
        return property(get, put)
    timeCourseSelections = selections = _selections('sel', roadrunner.roadrunner.RoadRunner.timeCourseSelections)
    steadyStateSelections = _selections('ssel', roadrunner.roadrunner.RoadRunner.steadyStateSelections)
    del _selections

    def _cma_set(self, v):
        if v and '_lz' in self.__dict__:
            self._load_rr()  # ponytail: libstruct's species order/moieties; Rust's own choice can differ
        roadrunner.roadrunner.RoadRunner.conservedMoietyAnalysis.fset(self, v)
    conservedMoietyAnalysis = property(roadrunner.roadrunner.RoadRunner.conservedMoietyAnalysis.fget, _cma_set)

    def getValue(self, *args):
        lz = self.__dict__.get('_lz')
        if lz is not None and len(args) == 1:
            v = self._lz_get(lz, args[0])
            if v is not None:
                return v
        return self._super().getValue(*args)
    getValue.__doc__ = roadrunner.roadrunner.RoadRunner.getValue.__doc__

    def setValue(self, *args):
        lz = self.__dict__.get('_lz')
        if lz is not None and len(args) == 2 and self._lz_set(lz, *args):
            return
        return self._super().setValue(*args)
    setValue.__doc__ = roadrunner.roadrunner.RoadRunner.setValue.__doc__

    def __getitem__(self, key):
        return self.getValue(key)

    def __setitem__(self, key, value):
        self.setValue(key, value)

    def __getattr__(self, name):  # lazy: model ids as attributes (roadrunner's dynamic properties)
        lz = self.__dict__.get('_lz')
        if lz is not None:
            if name in lz['attr']:
                return self.getValue(lz['attr'][name])
            raise AttributeError(name)
        return roadrunner.roadrunner.RoadRunner.__getattr__(self, name)

    def __setattr__(self, name, value):
        lz = self.__dict__.get('_lz')
        if lz is not None and name in lz['attr']:
            return self.setValue(lz['attr'][name], value)
        object.__setattr__(self, name, value)

    def getCurrentTime(self):
        lz = self.__dict__.get('_lz')
        return self._super().getCurrentTime() if lz is None else lz['t']
    getCurrentTime.__doc__ = roadrunner.roadrunner.RoadRunner.getCurrentTime.__doc__

    def isModelLoaded(self):
        return '_lz' in self.__dict__ or roadrunner.roadrunner.RoadRunner.isModelLoaded(self)

    @staticmethod
    def _lz_doc(lz):  # the SBML document as roadrunner keeps it
        import libsbml
        doc = libsbml.readSBMLFromString(lz['sbml'])
        for rx in doc.getModel().getListOfReactions():
            for l in (rx.getListOfReactants(), rx.getListOfProducts()):
                for sr in l:
                    if not sr.isSetStoichiometry() and not sr.isSetStoichiometryMath():
                        sr.setStoichiometry(sr.getStoichiometry())
        return doc

    def getSBML(self):
        lz = self.__dict__.get('_lz')
        if lz is None or lz['srefs']:
            return self._super().getSBML()
        import libsbml
        return libsbml.writeSBMLToString(self._lz_doc(lz))
    getSBML.__doc__ = roadrunner.roadrunner.RoadRunner.getSBML.__doc__

    def getCurrentSBML(self):
        lz = self.__dict__.get('_lz')
        if lz is None or lz['srefs']:
            return self._super().getCurrentSBML()
        import libsbml  # the SBML with the current values, as roadrunner writes it
        doc = self._lz_doc(lz)
        m, md = doc.getModel(), lz['model']
        v = dict(zip(md.getGlobalParameterIds(), md.getGlobalParameterValues()))
        v.update(zip(md.getCompartmentIds(), md.getCompartmentVolumes()))
        for p in m.getListOfParameters():
            p.setValue(v[p.getId()])
        for c in m.getListOfCompartments():
            c.setSize(v[c.getId()])
        for ids, conc, amt in ((md.getFloatingSpeciesIds(), md.getFloatingSpeciesConcentrations(), md.getFloatingSpeciesAmounts()),
                               (md.getBoundarySpeciesIds(), md.getBoundarySpeciesConcentrations(), md.getBoundarySpeciesAmounts())):
            for s, c, a in zip(ids, conc, amt):
                x = m.getSpecies(s)
                x.setInitialConcentration(c) if x.isSetInitialConcentration() else x.setInitialAmount(a)
        while m.getNumInitialAssignments():
            m.removeInitialAssignment(0)
        return libsbml.writeSBMLToString(doc)
    getCurrentSBML.__doc__ = roadrunner.roadrunner.RoadRunner.getCurrentSBML.__doc__

    def getFullStoichiometryMatrix(self):
        lz = self.__dict__.get('_lz')
        if lz is None or lz['srefs']:
            return self._super().getFullStoichiometryMatrix()
        N = np.array(lz['info']['N'], dtype=float)
        for i, j, sg, sid in lz['nref']:
            N[i, j] += sg * (self.getValue(sid) - lz['info']['values'][sid])
        return _named(N, lz['info']['rows'], lz['st'].rids)
    getFullStoichiometryMatrix.__doc__ = roadrunner.roadrunner.RoadRunner.getFullStoichiometryMatrix.__doc__

    def setIntegrator(self, name):
        if name not in ('cvode', 'gillespie'):  # other integrators need roadrunner's model
            self._load_rr()
        roadrunner.roadrunner.RoadRunner.setIntegrator(self, name)
    setIntegrator.__doc__ = roadrunner.roadrunner.RoadRunner.setIntegrator.__doc__

    integrator = property(roadrunner.roadrunner.RoadRunner.getIntegrator, setIntegrator)

    def getRatesOfChange(self):
        lz = self.__dict__.get('_lz')
        if lz is None or not lz['st'].plain:
            return self._super().getRatesOfChange()
        return np.array([self._lz_get(lz, s + "'") for s in lz['st'].ids])
    getRatesOfChange.__doc__ = roadrunner.roadrunner.RoadRunner.getRatesOfChange.__doc__

    def getRateRuleIds(self):
        lz = self.__dict__.get('_lz')
        return self._super().getRateRuleIds() if lz is None else list(lz['info']['rate_rule_ids'])
    getRateRuleIds.__doc__ = roadrunner.roadrunner.RoadRunner.getRateRuleIds.__doc__

    # ---------------------------------------------------------------------
    # Turbo (Rust engine) overrides. Each public method runs its Rust path through _run, which
    # records the outcome in lastRun and falls back to roadrunner when the path raises _Fallback.
    # ---------------------------------------------------------------------
    @property
    def lastRun(self):
        """What answered the last turbo-capable call: dict(method, engine='rust'|'roadrunner',
        reason=None or why roadrunner ran; for Rust time courses also solver, bdf_steps,
        rodas_steps, rejected, switched, atol). None before any such call."""
        r = self.__dict__.get('_last_run')
        return None if r is None else dict(r)

    def _run(self, method, rust, fallback):
        """rust() unless it raises _Fallback, then fallback(); records lastRun either way."""
        try:
            res = rust()
        except _Fallback as e:
            _no_fallback('%s: %s' % (method, e))
            self.__dict__['_last_run'] = dict(method=method, engine='roadrunner', reason=str(e))
            return fallback()
        info = self.__dict__.pop('_run_info', None) or {}
        self.__dict__['_last_run'] = dict(method=method, engine='rust', reason=None, **info)
        return res

    def _rs(self):
        """The _RsModel for the current model; raises _Fallback (turbo off/unsupported) with the reason."""
        lz = self.__dict__.get('_lz')
        if lz is not None:
            if _turbo:
                return lz['st']
            self._load_rr()
        if not _turbo:
            raise _Fallback('turbo is off' if _has_rs else 'roadrunner_rs is not installed')
        mdl = self.model
        if mdl is None:
            raise _Fallback('no model loaded')
        c = self.__dict__.get('_rs_cache')
        # rebuild on a CMA switch or when a compartment volume baked into the Rust model changes
        cma = self.conservedMoietyAnalysis
        if c is None or c[0] != (cma, mdl.getCompartmentVolumes()[c[2]].tobytes()):
            try:
                st = self._rs_build()
            except _Fallback as e:
                st = e
            cc = st.cconst if isinstance(st, _RsModel) else list(range(self.model.getNumCompartments()))
            c = self.__dict__['_rs_cache'] = ((cma, mdl.getCompartmentVolumes()[cc].tobytes()), st, cc)
        if isinstance(c[1], _Fallback):
            raise c[1]
        return c[1]

    def _rs_or_none(self):
        try:
            return self._rs()
        except _Fallback:
            return None

    def _rs_build(self):
        from ..analysis.rrs_bridge import from_roadrunner
        info = {}
        try:
            m = from_roadrunner(self, info)
        except (NotImplementedError, ValueError) as e:
            raise _Fallback('model not supported by the Rust engine: %s' % e)
        mdl = self.model
        ids = [list(f()) for f in (mdl.getFloatingSpeciesIds, mdl.getBoundarySpeciesIds, mdl.getGlobalParameterIds,
                                   mdl.getCompartmentIds, mdl.getReactionIds)]
        # steady state: roadrunner refuses events, and time-dependent rates have none
        nss = 'symbols/time' in self.getCurrentSBML() or mdl.getNumEvents() > 0
        return _RsModel(m, info, ids, mdl.getCompartmentVolumes(), nss)

    def _rs_x(self, st, mdl=None):
        """Current values of st.keys: parameters, compartments, boundary species (SBML units), species conc."""
        lz = self.__dict__.get('_lz')
        if lz is not None:
            return lz['x'].copy()
        mdl = self.model if mdl is None else mdl
        b = [mdl.getBoundarySpeciesConcentrations(), mdl.getBoundarySpeciesAmounts()] if st.bidx else []
        return np.concatenate([mdl.getGlobalParameterValues(), mdl.getCompartmentVolumes(), *b,
                               mdl.getFloatingSpeciesConcentrations()])[st.xidx]

    def _rs_sync(self, st, d, mdl=None, plan='sync'):
        """Write Rust values (dict id -> value; species as concentrations) back into the roadrunner model."""
        arr = isinstance(d, np.ndarray)  # all outputs, in output_ids order (else a dict id -> value)
        lz = self.__dict__.get('_lz')
        if lz is not None:  # no roadrunner model yet: the state is lz['x']
            if arr:
                lz['x'][:] = d[lz['xo']]
            else:
                for k, v in d.items():
                    if k in st.kpos:
                        lz['x'][st.kpos[k]] = v
            return
        mdl = self.model if mdl is None else mdl
        for f, idx, names, oidx in getattr(st, plan):
            if arr and f not in ('amount', 'conc'):
                getattr(mdl, f)(idx, d[oidx])
                continue
            if arr:
                have = list(zip(idx, names, d[oidx]))
            else:
                have = [(i, n, d[n]) for i, n in zip(idx, names) if n in d]
            if not have:
                continue
            if f in ('amount', 'conc'):
                for _, n, x in have:
                    roadrunner.RoadRunner.setValue(self, n if f == 'amount' else '[%s]' % n, float(x))
            else:
                i, _, x = zip(*have)
                getattr(mdl, f)(list(i), np.array(x, dtype=float))

    def _rs_tol(self, solver):
        """(rtol, atol) of a solver, or _Fallback for settings Rust doesn't take."""
        try:
            return float(solver.getValue('relative_tolerance')), float(solver.getValue('absolute_tolerance'))
        except TypeError:
            raise _Fallback('per-species tolerance vectors')

    def _rs_integrator(self, integ):
        """The integrator's name if Rust honours it with its current settings, else _Fallback."""
        name = integ.getName()
        if name not in _RS_INTEGRATORS:
            raise _Fallback('integrator %s' % name)
        keys, vals = _RS_INTEGRATORS[name]
        cfg = integ.getSettingsMap()
        bad = [k for k, v in zip(keys, vals) if cfg.get(k) != v]
        if bad:
            raise _Fallback('%s settings %s' % (name, ', '.join(bad)))
        return name

    def _rs_advance(self, st, times, ores, integ):
        """Run the Rust model from the current state over times and write the final state back
        (species, event-changed values, time). Returns the outputs at times."""
        rtol, atol = self._rs_tol(integ)
        mdl = self.model
        Y, final = st.advance(self._rs_x(st, mdl).tolist(), times, ores, rtol, atol,
                              st.atol_scale(mdl.getFloatingSpeciesInitConcentrations()))
        if not np.isfinite(Y).all() and (~np.isfinite(Y) & np.isfinite(Y[:1])).any():  # went non-finite during the run
            raise _Fallback('non-finite Rust output')
        self._rs_sync(st, final, mdl)
        mdl.setTime(times[-1])
        self.__dict__['_run_info'] = st.stats
        return Y

    def _rs_mcast(self):
        st = self._rs()
        if not st.mca:
            raise _Fallback('rate rules, events or time-dependent rates (no steady state in roadrunner either)'
                            if st.plain else 'rate-rule or non-plain species')
        return st

    def simulate(self, start=None, end=None, points=None, selections=None, output_file=None, steps=None, times=None):
        # Turbo fast path for simulate(start, end, points or steps[, selections]) with CVODE: what
        # roadrunner's wrapper does (its simulate options are private and reset by every call) without
        # its per-call overhead. Anything else, or a model/setting Rust can't run, takes the wrapper
        # (which also turns gillespie's variable_step_size off for such calls, and calls _simulate).
        def wrapper():
            return super(ExtendedRoadRunner, self).simulate(start, end, points, selections, output_file, steps, times)
        if not (_turbo and output_file is None and times is None and start is not None and end is not None
                and (points is None) != (steps is None) and type(points) is not list):
            return wrapper()
        n = points - 1 if steps is None else steps
        integ = self.integrator
        sel = None
        try:
            if n < 1 or integ.getName() != 'cvode':
                raise _Fallback('wrapper')
            st = self._rs()
            sel = list(self.timeCourseSelections)
        except _Fallback:
            return wrapper()
        # new selections go through roadrunner's wrapper once (it validates and sets them)
        if not (selections is None or (isinstance(selections, (list, tuple)) and sel == list(selections))):
            return wrapper()

        def rust():
            res = self.__dict__['_rs_data'] = self._rs_simulate(st, _Opts(start=start, end=end, steps=n), sel, integ)
            return res
        return self._run('simulate', rust, wrapper)
    simulate.__doc__ = roadrunner.RoadRunner.simulate.__doc__

    def _simulate(self, o):
        def rust():
            return self._rs_simulate(self._rs(), o)

        def fallback():
            self.__dict__['_rs_data'] = None
            st = self._rs_or_none()
            if st is not None:
                st.reset_events()  # roadrunner moves the event triggers itself; Rust restarts them next time
            return self._super()._simulate(o)
        res = self._run('simulate', rust, fallback)
        if self.__dict__.get('_last_run', {}).get('engine') == 'rust':
            self.__dict__['_rs_data'] = res
        return res

    def _rs_simulate(self, st, o, sel=None, integ=None):
        integ = self.integrator if integ is None else integ
        if o.output_file:
            raise _Fallback('output_file')
        name = self._rs_integrator(integ)
        times = np.asarray(o.times, dtype=float) if len(o.times) else st.times(o.start, o.end, o.steps)
        sel = list(self.timeCourseSelections) if sel is None else sel
        cols = st.cols(sel)
        if cols is None:
            raise _Fallback('selections %s' % sel)
        if o.reset_model:
            self.reset()
        pre = cols.init_values(self)  # init(...) values before the run
        if name == 'gillespie':
            Y = self._rs_ssa(st, o, times, cols.outs)
        else:
            Y = self._rs_advance(st, times, cols.ores, integ)
        return _named(cols.apply(Y, pre), cols=sel)

    def _rs_ssa(self, st, o, times, outs):
        """Gillespie in Rust: fixed output points, or every reaction firing when variable_step_size.
        Writes the final state back; returns the output rows."""
        integ = self.integrator
        # roadrunner draws successive runs from one seeded stream: same seed -> same sequence of runs
        seed = int(integ.getValue('seed')) % 2 ** 64
        if seed == 2 ** 64 - 1 and '_lz' in self.__dict__:  # unset: roadrunner seeds from the clock at load
            seed = self.__dict__.setdefault('_rs_auto_seed', np.random.SeedSequence().entropy % (2 ** 64 - 1))
        c = self.__dict__.get('_rs_seed')
        if c is None or c[0] != seed:
            c = self.__dict__['_rs_seed'] = [seed, 0]
        c[1] += 1
        run_seed = (seed + c[1] * 0x9E3779B97F4A7C15) % 2 ** 64
        m, x, k = st.m, self._rs_x(st), st.keys
        try:
            if integ.getValue('variable_step_size'):
                R = roadrunner_rs.gillespie_steps(m, k, x.tolist(), o.start, o.end, run_seed, outputs=outs + k,
                                                  max_rows=int(integ.getValue('max_output_rows')))
                Y, tend = R[:, 1:], R[-1, 0]
            else:
                Y, tend = m.gillespie_batch(k, x[None], times, seed=run_seed, threads=1, outputs=outs + k)[0], times[-1]
        except ValueError:
            raise _Fallback('rate rules in a stochastic simulation')
        if not np.isfinite(Y).all() and (~np.isfinite(Y) & np.isfinite(Y[:1])).any():  # went non-finite during the run
            raise _Fallback('non-finite Rust output')
        st.reset_events()  # the next ODE run starts its events fresh
        self._rs_sync(st, dict(zip(k, Y[-1, len(outs):])))
        self.model.setTime(tend)
        self.__dict__['_run_info'] = dict(solver='gillespie')
        return Y[:, :len(outs)]

    def timeSeriesSensitivities(self, start, stop, num, params=(), species=(), k=0):
        return self._run('timeSeriesSensitivities',
                         lambda: self._rs_tss(self._rs(), start, stop, num, list(params), list(species), k),
                         lambda: self._super().timeSeriesSensitivities(start, stop, num, list(params), list(species), k))
    timeSeriesSensitivities.__doc__ = roadrunner.RoadRunner.timeSeriesSensitivities.__doc__

    def _rs_tss(self, st, start, stop, num, params, species, k):
        """Forward sensitivities in Rust: exact symbolic derivatives, parameters integrated in parallel."""
        if k not in (0, 1):
            raise _Fallback('sensitivity order k=%s' % k)
        cma = self.conservedMoietyAnalysis
        cs = list(self.getConservedMoietyIds()) if cma else []
        # roadrunner: rows = global parameters (+ moiety totals, + rate-rule targets as zero rows);
        # columns = rate-rule states, then (independent) floating species, in model order
        rr_ids = [v for v in self.getRateRuleIds() if v in st.states]
        base = rr_ids + (list(self.getIndependentFloatingSpeciesIds()) if cma else st.ids)
        pkeys = set(st.keys[:st.off]) - set(st.bids)
        if any(p not in pkeys and p not in cs for p in params) or \
                any(s in st.general or s not in st.states for s in base):
            raise _Fallback('sensitivity parameters/species the Rust model does not hold')
        params, species = params or st.gids, [s for s in base if s in species or not species]
        if '_lz' in self.__dict__:  # roadrunner makes its sensitivity solver with the model: still the defaults
            rtol, atol = 1e-6, 1e-12
        else:
            rtol, atol = self._rs_tol(self.getSensitivitySolver())
        zero = set(rr_ids)
        rp = [p for p in params if p not in cs and p not in zero]
        dep = list(self.getDependentFloatingSpeciesIds()) if any(p in cs for p in params) else []
        times = np.linspace(start, stop, int(num))
        try:
            X, S = roadrunner_rs.time_sensitivities(st.m, st.keys, self._rs_x(st).tolist(), times, rp + dep,
                                                    rtol=rtol, atol=atol, order=k)
        except ValueError:
            raise _Fallback('Rust sensitivity integration failed')
        rows = []
        if dep:  # moiety total T_k moves the dependent species: d conc_dep / d T = inv(Gamma[:, dep]) / V_dep
            G = self.getConservationMatrix()
            gcols = list(G.colnames)
            C = np.linalg.inv(np.asarray(G)[:, [gcols.index(d) for d in dep]]) / st.V[[st.ids.index(d) for d in dep]][:, None]
        for p in params:
            if p in zero:
                rows.append(np.zeros((len(times), S.shape[2])))
            elif p in cs:
                rows.append(np.einsum('tdx,d->tx', S[:, len(rp):], C[:, cs.index(p)]))
            else:
                rows.append(S[:, rp.index(p)])
        S = np.stack(rows, axis=1) if rows else S[:, :0]
        if k:
            S[0] = 0.0  # roadrunner reports a zero first row for k = 1 too
        self._rs_sync(st, dict(zip(st.states, X[-1])))
        self.model.setTime(times[-1])
        idx = [st.states.index(s) for s in species]
        V = np.array([st.V[st.ids.index(s)] if s in st.ids else 1.0 for s in species])
        return times, S[:, :, idx] * V, params, species  # roadrunner reports amount sensitivities

    def oneStep(self, currentTime, stepSize, reset=True):
        def rust():
            st, integ = self._rs(), self.integrator
            if self._rs_integrator(integ) != 'cvode':
                raise _Fallback('integrator %s' % integ.getName())
            t = np.array([currentTime, currentTime + stepSize])
            self._rs_advance(st, t, None, integ)
            return currentTime + stepSize
        return self._run('oneStep', rust, lambda: self._super().oneStep(currentTime, stepSize, reset))

    def simulateBatch(self, params, start, end, points, selections=None, threads=None):
        """Many time courses at once on all CPU cores (Rust engine), one per parameter set.

        :param params: dict ``{'k1': [...], 'k2': [...]}`` of equal-length sequences (zipped, not a
            grid), or a pandas DataFrame. Ids can be parameters, compartments, boundary or floating
            species (concentrations), as in ``r[id] = value``.
        :param start, end, points: as in ``simulate``.
        :param selections: output columns (default: timeCourseSelections).
        :param threads: CPU threads (default: all cores).
        :returns: array (n sets, points, len(selections)). Every run starts from the current state;
            the model itself is not changed. Runs Rust cannot do (or that fail) run in roadrunner.
        """
        if hasattr(params, 'to_dict'):
            params = params.to_dict('list')
        names = list(params)
        vals = np.array([np.atleast_1d(np.asarray(params[k], dtype=float)) for k in names]).T.reshape(-1, len(names))
        sel = list(self.timeCourseSelections if selections is None else selections)
        times = np.linspace(start, end, points)

        def rust():
            st, integ = self._rs(), self.integrator
            if self._rs_integrator(integ) != 'cvode':
                raise _Fallback('integrator %s' % integ.getName())
            cols = st.cols(sel)
            if cols is None or cols.init:
                raise _Fallback('selections %s' % sel)
            bad = [k for k in names if k.strip('[]') not in st.kpos]
            if bad:
                raise _Fallback('not a Rust parameter or state: %s' % ', '.join(bad))
            rtol, atol = self._rs_tol(integ)
            mdl = self.model
            X = np.repeat(self._rs_x(st, mdl)[None], len(vals), axis=0)
            for j, k in enumerate(names):  # species ids are concentrations, as r['[S]'] / r.S for conc
                X[:, st.kpos[k.strip('[]')]] = vals[:, j]
            Y = st.m.simulate_batch(st.keys, X, times, threads=threads, rtol=rtol, atol=atol, outputs=cols.outs,
                                    atol_scale=st.atol_scale(mdl.getFloatingSpeciesInitConcentrations()))
            out = np.stack([cols.apply(y) for y in Y]) if len(Y) else np.zeros((0, points, len(sel)))
            for i in np.flatnonzero(~np.isfinite(out).all(axis=(1, 2))):  # failed runs: roadrunner
                out[i] = self._rr_one(names, vals[i], start, end, points, sel)
            return out

        def fallback():
            return np.stack([self._rr_one(names, v, start, end, points, sel) for v in vals]) if len(vals) \
                else np.zeros((0, points, len(sel)))
        return self._run('simulateBatch', rust, fallback)

    def _rr_one(self, names, v, start, end, points, sel):
        """One simulateBatch run in plain roadrunner, on a copy of the current model (this one unchanged)."""
        r = roadrunner.roadrunner.RoadRunner(self.getCurrentSBML())
        r.setIntegrator(self.integrator.getName())
        for k, x in self.integrator.getSettingsMap().items():
            try:
                r.integrator.setValue(k, x)
            except Exception:  # read-only or engine-specific setting
                pass
        sp = set(r.model.getFloatingSpeciesIds()) | set(r.model.getBoundarySpeciesIds())
        for k, x in zip(names, v):
            r.setValue('[%s]' % k if k in sp else k, float(x))  # species ids are concentrations, as in Rust
        return np.asarray(r.simulate(start, end, points, sel))

    def getSimulationData(self):
        d = self.__dict__.get('_rs_data')
        return self._super().getSimulationData() if d is None else d

    def steadyState(self, dict=None):
        def rust():
            st = self._rs()
            if dict is not None:
                raise _Fallback('steady-state options dict')
            if st.nss:
                raise _Fallback('events or time-dependent rates')
            if self.steadyStateSolver.getName() != 'nleq2':
                raise _Fallback('steady-state solver %s' % self.steadyStateSolver.getName())
            S, v, _, ok = st.m.steady_state_batch(st.keys, self._rs_x(st)[None])
            if not ok[0]:
                raise _Fallback('Rust Newton found no steady state')
            self._rs_sync(st, {k: val for k, val in zip(st.states, S[0])})
            dx = (st.N @ v[0][st.rx] / st.V) if st.plain else 0.0
            return float(np.linalg.norm(dx))
        return self._run('steadyState', rust, lambda: self._super().steadyState(dict))

    def _rs_mca(self):
        """(st, unscaled CCC, unscaled FCC, v, S) at steady state, model left there (as roadrunner does)."""
        st = self._rs_mcast()
        res = st.m.mca_batch(st.keys, self._rs_x(st)[None])
        if not res['ok'][0]:
            raise _Fallback('Rust found no steady state for MCA')
        sp, rx = st.sp, st.rx
        self._rs_sync(st, dict(zip(st.states, res['S'][0])))
        return st, res['uccc'][0][sp][:, rx], res['ufcc'][0][rx][:, rx], res['v'][0][rx], res['S'][0][sp]

    def getUnscaledConcentrationControlCoefficientMatrix(self):
        def rust():
            st, C, _, _, _ = self._rs_mca()
            return _named(C, st.ids, st.rids)
        return self._run('getUnscaledConcentrationControlCoefficientMatrix', rust,
                         lambda: self._super().getUnscaledConcentrationControlCoefficientMatrix())

    def getScaledConcentrationControlCoefficientMatrix(self):
        def rust():
            st, C, _, v, S = self._rs_mca()
            return _named(C * v[None, :] / S[:, None], st.ids, st.rids)
        return self._run('getScaledConcentrationControlCoefficientMatrix', rust,
                         lambda: self._super().getScaledConcentrationControlCoefficientMatrix())

    def getUnscaledFluxControlCoefficientMatrix(self):
        def rust():
            st, _, F, _, _ = self._rs_mca()
            return _named(F, st.rids, st.rids)
        return self._run('getUnscaledFluxControlCoefficientMatrix', rust,
                         lambda: self._super().getUnscaledFluxControlCoefficientMatrix())

    def getScaledFluxControlCoefficientMatrix(self):
        def rust():
            st, _, F, v, _ = self._rs_mca()
            vi = np.where(np.abs(v) < self.getFluxThreshold(), 0.0, v)  # roadrunner zeroes near-zero flux rows
            with np.errstate(divide='ignore', invalid='ignore'):
                out = np.where(vi[:, None] != 0, F * v[None, :] / vi[:, None], 0.0)
            return _named(out, st.rids, st.rids)
        return self._run('getScaledFluxControlCoefficientMatrix', rust,
                         lambda: self._super().getScaledFluxControlCoefficientMatrix())

    def _rs_rd(self):
        """(st, v, dv/dS, S) for the reactions at the current state (exact derivatives)."""
        st = self._rs()
        if not st.plain:
            raise _Fallback('rate-rule or non-plain species')
        x = self._rs_x(st)
        if st.rd is None:  # compiled once per Rust model
            st.rd = roadrunner_rs.RateDerivatives(st.m)
        v, dx, _, _ = st.rd(st.keys, x.tolist(), self.model.getTime())
        return st, v[st.rx], dx[st.rx][:, st.sp], x[st.off:]

    def getUnscaledElasticityMatrix(self):
        def rust():
            st, _, E, _ = self._rs_rd()
            return _named(E, st.rids, st.ids)
        return self._run('getUnscaledElasticityMatrix', rust, lambda: self._super().getUnscaledElasticityMatrix())

    def getScaledElasticityMatrix(self):
        def rust():
            st, v, E, S = self._rs_rd()
            with np.errstate(divide='ignore', invalid='ignore'):
                return _named(E * S[None, :] / v[:, None], st.rids, st.ids)
        return self._run('getScaledElasticityMatrix', rust, lambda: self._super().getScaledElasticityMatrix())

    def getSteadyStateValues(self):
        if self._rs_or_none() is None:
            return self._super().getSteadyStateValues()
        self.steadyState()  # Rust when it can, roadrunner otherwise
        return np.array([self.getValue(s) for s in self.steadyStateSelections], dtype=float)

    def getSteadyStateValuesNamedArray(self):
        if self._rs_or_none() is None:
            return self._super().getSteadyStateValuesNamedArray()
        return _named(self.getSteadyStateValues()[None], cols=self.steadyStateSelections)

    def steadyStateNamedArray(self, dict=None):
        if self._rs_or_none() is None or dict is not None:
            return self._super().steadyStateNamedArray(dict)
        return self.getSteadyStateValuesNamedArray()

    def mcaSteadyState(self, *args):
        return self._run('mcaSteadyState', lambda: (self._rs_mcast(), self.steadyState())[1],
                         lambda: self._super().mcaSteadyState(*args))

    def _rs_jac(self):
        """(st, Jacobian of the floating species in concentrations) at the current state and time."""
        st = self._rs_mcast()
        sp = st.sp
        J = np.asarray(st.m.jacobian(st.keys, self._rs_x(st).tolist(), self.model.getTime()))
        return st, J[sp][:, sp]

    def getFullJacobian(self):
        def rust():
            st, J = self._rs_jac()
            return _named(J, st.ids, st.ids)
        return self._run('getFullJacobian', rust, lambda: self._super().getFullJacobian())

    def _rs_rjac(self):
        st, J = self._rs_jac()
        if not self.conservedMoietyAnalysis:
            return st.ids, J
        if st.ind is None:  # roadrunner's (libstruct's) choice of independent species
            Lm = self.getLinkMatrix()
            ok = set(st.ids) <= set(Lm.rownames)
            st.L = np.array(Lm)[[list(Lm.rownames).index(s) for s in st.ids]] if ok else None
            st.ind = list(Lm.colnames)
        if st.L is None:
            raise _Fallback("link matrix rows don't cover the floating species")
        # d(dependent conc)/d(independent conc): the link matrix is in amounts
        V, ind = st.V, [st.ids.index(s) for s in st.ind]
        Lc = st.L / V[:, None] * V[ind][None, :]
        return st.ind, (J @ Lc)[ind]

    def getReducedJacobian(self, h=-1.0):
        def rust():
            ids, J = self._rs_rjac()
            return _named(J, ids, ids)
        return self._run('getReducedJacobian', rust, lambda: self._super().getReducedJacobian(h))

    def getFullEigenValues(self):
        return self._run('getFullEigenValues', lambda: np.linalg.eigvals(self._rs_jac()[1]).astype(complex),
                         lambda: self._super().getFullEigenValues())

    def getReducedEigenValues(self):
        return self._run('getReducedEigenValues', lambda: np.linalg.eigvals(self._rs_rjac()[1]).astype(complex),
                         lambda: self._super().getReducedEigenValues())

    def getFullEigenValuesNamedArray(self):
        def rust():
            st, J = self._rs_jac()
            e = np.linalg.eigvals(J)
            return _named(np.column_stack([e.real, e.imag]), st.ids, ['real', 'imaginary'])
        return self._run('getFullEigenValuesNamedArray', rust, lambda: self._super().getFullEigenValuesNamedArray())

    def getReducedEigenValuesNamedArray(self):
        def rust():
            ids, J = self._rs_rjac()
            e = np.linalg.eigvals(J)
            return _named(np.column_stack([e.real, e.imag]), ids, ['real', 'imaginary'])
        return self._run('getReducedEigenValuesNamedArray', rust, lambda: self._super().getReducedEigenValuesNamedArray())

    def _rs_ucc(self, variableName, parameterName):
        """(unscaled CC, parameter value, variable value now, is flux) with roadrunner's 5-point stencil."""
        st = self._rs_mcast()
        var = variableName.replace('[', '').replace(']', '')
        if st.kpos.get(parameterName, st.off) >= st.off or (var not in st.ids and var not in st.rids):
            raise _Fallback('species / moiety totals as parameters')
        x = self._rs_x(st)
        j = st.kpos[parameterName]
        p0 = x[j]
        h = self.getDiffStepSize() * p0
        if abs(h) < 1e-12:
            h = self.getDiffStepSize()
        X = np.repeat(x[None], 4, axis=0)
        X[:, j] += [h, 2 * h, -h, -2 * h]
        S, v, _, ok = st.m.steady_state_batch(st.keys, X)
        if not ok.all():
            raise _Fallback('Rust found no steady state for a stencil point')
        if var in st.ids:
            i = st.ids.index(var)
            f, now = S[:, st.sp[i]], self.model.getFloatingSpeciesConcentrations()[i]
        else:
            i = st.rids.index(var)
            f, now = v[:, st.rx[i]], self.model.getReactionRates()[i]
        return (f[3] + 8 * f[0] - 8 * f[2] - f[1]) / (12 * h), p0, now, var in st.rids

    def getuCC(self, variableName, parameterName):
        return self._run('getuCC', lambda: self._rs_ucc(variableName, parameterName)[0],
                         lambda: self._super().getuCC(variableName, parameterName))

    def getCC(self, variableName, parameterName):
        def rust():
            u, p, now, flux = self._rs_ucc(variableName, parameterName)
            if flux and abs(now) < self.getFluxThreshold():
                return 0.0
            with np.errstate(divide='ignore', invalid='ignore'):
                return float(np.float64(u) * p / now)
        return self._run('getCC', rust, lambda: self._super().getCC(variableName, parameterName))

    # ---------------------------------------------------------------------
    # Model access
    # ---------------------------------------------------------------------
    # def getBoundarySpeciesConcentrations(self):
    #     return self.model.getBoundarySpeciesConcentrations()
    # getBoundarySpeciesConcentrations.__doc__ = roadrunner.ExecutableModel.getBoundarySpeciesConcentrations.__doc__

    # These model functions are attached after class creation
    _model_functions = [
        'getBoundarySpeciesIds',
        'getBoundarySpeciesConcentrationIds',
        'getBoundarySpeciesConcentrations',
        'getBoundarySpeciesAmounts',
        'getNumBoundarySpecies',

        'getFloatingSpeciesIds',
        'getFloatingSpeciesConcentrationIds',
        'getFloatingSpeciesConcentrations',
        'getFloatingSpeciesAmounts',
        'getNumFloatingSpecies',

        'getGlobalParameterIds',
        'getGlobalParameterValues',
        'getNumGlobalParameters',

        'getCompartmentIds',
        'getCompartmentVolumes',
        'getNumCompartments',

        'getConservedMoietyValues',
        'getNumConservedMoieties',

        'getReactionIds',
        'getReactionRates',
        'getNumReactions',
        'getNumEvents',
        'getNumRateRules'
     ]

    # ---------------------------------------------------------------------
    # Jarnac compatibility layer
    # ---------------------------------------------------------------------
    def fjac(self):
        return self.getFullJacobian()
    fjac.__doc__ = roadrunner.RoadRunner.getFullJacobian.__doc__

    def sm(self):
        return self.getFullStoichiometryMatrix()
    sm.__doc__ = roadrunner.RoadRunner.getFullStoichiometryMatrix.__doc__

    def rs(self):
        return self.model.getReactionIds()
    rs.__doc__ = roadrunner.ExecutableModel.getReactionIds.__doc__

    def fs(self):
        return self.model.getFloatingSpeciesIds()
    fs.__doc__ = roadrunner.ExecutableModel.getFloatingSpeciesIds.__doc__

    def bs(self):
        return self.model.getBoundarySpeciesIds()
    bs.__doc__ = roadrunner.ExecutableModel.getBoundarySpeciesIds.__doc__

    def ps(self):
        return self.model.getGlobalParameterIds()
    ps.__doc__ = roadrunner.ExecutableModel.getGlobalParameterIds.__doc__

    def vs(self):
        return self.model.getCompartmentIds()
    vs.__doc__ = roadrunner.ExecutableModel.getCompartmentIds.__doc__

    def dv(self):
        return self.getRatesOfChange()
    dv.__doc__ = roadrunner.RoadRunner.getRatesOfChange.__doc__

    def rv(self):
        return self.model.getReactionRates()
    rv.__doc__ = roadrunner.ExecutableModel.getReactionRates.__doc__

    def sv(self):
        return self.model.getFloatingSpeciesConcentrations()
    sv.__doc__ = roadrunner.ExecutableModel.getFloatingSpeciesConcentrations.__doc__

    # ---------------------------------------------------------------------
    # Export Utilities
    # ---------------------------------------------------------------------
    def __getSBML(self, current):
        if current is True:
            return self.getCurrentSBML()
        else:
            return self.getSBML()

    def getAntimony(self, current=False, removeFunctionDefinitions=None):
        """ Antimony string of the original model loaded into roadrunner.

        :param current: return current model state
        :param removeFunctionDefinitions: remove function definitions
        :type current: bool
        :return: Antimony
        :rtype: str
        """
        sbml = self.__getSBML(current)
        from .. import sbmlToAntimony
        return sbmlToAntimony(sbml, removeFunctionDefinitions)

    def getCurrentAntimony(self, removeFunctionDefinitions=None):
        """ Antimony string of the current model state.

        See also: :func:`getAntimony`
        :returns: Antimony string
        :rtype: str
        """
        return self.getAntimony(current=True)

    def getCellML(self, current=False):
        """ CellML string of the original model loaded into roadrunner.

        :param current: return current model state
        :type current: bool
        :returns: CellML string
        :rtype: str
        """
        sbml = self.__getSBML(current)
        from .. import sbmlToCellML
        return sbmlToCellML(sbml)

    def getCurrentCellML(self):
        """ CellML string of current model state.

        See also: :func:`getCellML`
        :returns: CellML string
        :rtype: str
        """
        return self.getCellML(current=True)

    def getMatlab(self, current=False):
        """ Matlab string of the original model loaded into roadrunner.

        See also: :func:`getCurrentMatlab`
        :returns: Matlab string
        :rtype: str
        """
        try:
            from sbml2matlab import sbml2matlab
            sbml = self.__getSBML(current)
            return sbml2matlab.sbml2matlab(sbml)
        except ImportError:
            warnings.warn("'sbml2matlab' could not be imported, no support for Matlab code generation",
                          RuntimeWarning, stacklevel=2)
            return ""

    def getCurrentMatlab(self):
        """ Matlab string of current model state.

        :param current: return current model state
        :type current: bool
        :returns: Matlab string
        :rtype: str
        """
        return self.getMatlab(current=True)

    def exportToSBML(self, filePath, current=True):
        """ Save current model as SBML file.

        :param current: export current model state
        :type current: bool
        :param filePath: file path of SBML file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.__getSBML(current))

    def exportToAntimony(self, filePath, current=True):
        """ Save current model as Antimony file.

        :param current: export current model state
        :type current: bool
        :param filePath: file path of Antimony file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.getAntimony(current))

    def exportToCellML(self, filePath, current=True):
        """ Save current model as CellML file.

        :param current: export current model state
        :type current: bool
        :param filePath: file path of CellML file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.getCellML(current))

    def exportToMatlab(self, filePath, current=True):
        """ Save current model as Matlab file.
        To save the original model loaded into roadrunner use
        current=False.

        :param self: RoadRunner instance
        :type self: RoadRunner.roadrunner
        :param filePath: file path of Matlab file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.getMatlab(current))

    # ---------------------------------------------------------------------
    # Plotting Utilities
    # ---------------------------------------------------------------------
    def draw(self, **kwargs):
        """ Draws an SBMLDiagram of the current model.

        To set the width of the output plot provide the 'width' argument.
        Species are drawn as white circles (boundary species
        shaded in blue), reactions as grey squares.
        Currently only the drawing of medium-size networks is supported.
        """
        import shutil
        # PATH variable may not be present - cannot use os.environ['PATH']
        # but can use shutil.which for Python 3.3+
        if hasattr(shutil, 'which'):
            if shutil.which('dot') is None:
                warnings.warn("Graphviz is not installed in your machine or could not be found. 'draw' command cannot produce a diagram.",
                    Warning, stacklevel=2)
                return
        elif not 'PATH' in os.environ or any([ os.access( os.path.join( p, 'dot' ), os.X_OK ) for p in os.environ['PATH'].split( os.pathsep )]):
            warnings.warn("Graphviz is not installed in your machine or could not be found. 'draw' command cannot produce a diagram.",
                Warning, stacklevel=2)
            return

        from tellurium import SBMLDiagram
        diagram = SBMLDiagram(self.getSBML())
        diagram.draw(**kwargs)


    def plot(self, result=None, show=True,
             xlabel=None, ylabel=None, title=None, linewidth=2, xlim=None, ylim=None, logx=False, logy=False,
             xscale='linear', yscale='linear', grid=False, ordinates=None, tag=None, labels=None, 
             figsize=(6,4), savefig=None, dpi=80, alpha=1.0, **kwargs):
        """ Plot roadrunner simulation data.

        Plot is called with simulation data to plot as the first argument. If no data is provided the data currently
        held by roadrunner generated in the last simulation is used. The first column is considered the x axis and
        all remaining columns the y axis.
        If the result array has no names, then the current r.selections are used for naming. In this case the
        dimension of the r.selections has to be the same like the number of columns of the result array.

        Curves are plotted in order of selection (columns in result).

        In addition to the listed keywords plot supports all matplotlib.pyplot.plot keyword arguments,
        like color, alpha, linewidth, linestyle, marker, ...
        ::

            sbml = te.getTestModel('feedback.xml')
            r = te.loadSBMLModel(sbml)
            s = r.simulate(0, 100, 201)
            r.plot(s, loc="upper right", linewidth=2.0, lineStyle='-', marker='o', markersize=2.0,
                   alpha=0.8, title="Feedback Oscillation", xlabel="time", ylabel="concentration",
                   xlim=[0,100], ylim=[-1, 4])

        :param result: results data to plot (numpy array)
        :param show: show=True (default) shows the plot, use show=False to plot multiple simulations in one plot
        :type show: bool
        :param xlabel: x-axis label
        :type xlabel: str
        :param ylabel: y-axis label
        :type ylabel: str
        :param title: plot title
        :type title: str
        :param linewidth: linewidth of the plot
        :type linewidth: float
        :param xlim: limits on x-axis (tuple [start, end])
        :param ylim: limits on y-axis
        :param logx: use log scale for x-axis
        :type logx: bool
        :param logy: use log scale for y-axis
        :type logy: bool
        :param xscale: 'linear' or 'log' scale for x-axis
        :param yscale: 'linear' or 'log' scale for y-axis
        :param grid: show grid
        :type grid: bool
        :param ordinates: If supplied, only these selections will be plotted (see RoadRunner selections)
        :param tag: If supplied, all traces with the same tag will be plotted with the same color/style
        :param labels: 'id' to use species IDs
        :param figsize: If supplied, customize the size of the figure (width,height)
        :param savefig: If supplied, saves the figure to specified location
        :type savefig: str
        :param dpi: Change the dpi of the saved figure
        :type dpi: int
        :param alpha: Change the alpha value of the figure
        :type alpha: float
        :param kwargs: additional matplotlib keywords like marker, lineStyle, ...
        """
        if result is None:
            simData = self.getSimulationData()
            result = copy.copy(simData)
            result.colnames = simData.colnames

        from .. import getPlottingEngine

        if xlabel:
            kwargs['xlabel'] = xlabel
        if ylabel:
            kwargs['ylabel'] = ylabel
        if title:
            kwargs['title'] = title
        if linewidth:
            kwargs['linewidth'] = linewidth
        if xlim:
            kwargs['xlim'] = xlim
        if ylim:
            kwargs['ylim'] = ylim
        if logx:
            kwargs['logx'] = logx
        if logy:
            kwargs['logy'] = logy
        if xscale:
            kwargs['xscale'] = xscale
        if yscale:
            kwargs['yscale'] = yscale
        if grid:
            kwargs['grid'] = grid
        if ordinates:
            kwargs['ordinates'] = ordinates
        if tag:
            kwargs['tag'] = tag
        if labels:
            kwargs['labels'] = labels
        if figsize:
            kwargs['figsize'] = figsize
        if savefig:
            kwargs['savefig'] = savefig
        if dpi:
            kwargs['dpi'] = dpi
        if alpha:
            kwargs['alpha'] = alpha            

        # FIXME: provide the additional parameters to the plotting engine

        engine = getPlottingEngine()
        if show:
            # show the plot immediately
            engine.plotTimecourse(result, **kwargs)
        else:
            # otherwise, accumulate the traces
            engine.accumulateTimecourse(result, **kwargs)

    def show(self, reset=True):
        """ Show plot. """
        from .. import getPlottingEngine
        getPlottingEngine().show(reset=reset)

    # ---------------------------------------------------------------------
    # Stochastic Simulation Methods
    # ---------------------------------------------------------------------

    @property
    def seed(self):
        """ Getter for global Configuration seed. """
        return self.getSeed()

    @seed.setter
    def seed(self, value):
        """ Setter for model and global Configuration option seed. """
        return self.setSeed(value)

    def gillespie(self, *args, **kwargs):
        """ Run a Gillespie stochastic simulation.

        Sets the integrator to gillespie and performs simulation.
        ::

            rr = te.loada ('S1 -> S2; k1*S1; k1 = 0.1; S1 = 40')
            # Simulate from time zero to 40 time units using variable step sizes (classic Gillespie)
            result = rr.gillespie (0, 40)
            # Simulate on a grid with 10 points from start 0 to end time 40
            rr.reset()
            result = rr.gillespie (0, 40, 10)
            # Simulate from time zero to 40 time units using variable step sizes with given selection list
            # This means that the first column will be time and the second column species S1
            rr.reset()
            result = rr.gillespie (0, 40, selections=['time', 'S1'])
            # Simulate on a grid with 20 points from time zero to 40 time units
            # using the given selection list
            rr.reset()
            result = rr.gillespie (0, 40, 20, ['time', 'S1'])
            rr.plot(result)

        :param seed: seed for gillespie
        :type seed: int
        :param args: parameters for simulate
        :param kwargs: parameters for simulate
        :returns: simulation results
        """

        integratorName = self.integrator.getName()
        vss = self.integrator.variable_step_size
        if (integratorName != 'gillespie'):
            self.setIntegrator('gillespie')
            self.integrator.variable_step_size = True        
        if (len(args) > 3):
            self.integrator.variable_step_size = False
        if (len(args) == 3):
            if (type(args[2]) != list):
                self.integrator.variable_step_size = False
        elif ('points' in kwargs or 'steps' in kwargs):
            self.integrator.variable_step_size = False
        if "seed" in kwargs:
            self.integrator.setValue("seed", kwargs["seed"])
            del kwargs["seed"]

        s = self.simulate(*args, **kwargs)
        self.setIntegrator(integratorName)
        self.integrator.variable_step_size = vss
        return s


def _resets(f, everything):  # reset() / resetAll() / resetToOrigin() restart event triggers, as Rust does
    def g(self, *args, **kwargs):
        lz = self.__dict__.get('_lz')
        if lz is not None and not args and not kwargs and (everything or lz['reset_ok']):
            return self._lz_reset(lz, everything)  # lazily nothing can change the SBML: origin = all
        self._load_rr()
        if f.__name__ == 'resetToOrigin':
            self.__dict__.pop('_rs_cache', None)
        c = self.__dict__.get('_rs_cache')
        if c is not None and isinstance(c[1], _RsModel):
            c[1].reset_events()
        return f(self, *args, **kwargs)
    g.__name__, g.__doc__, g._route = f.__name__, f.__doc__, 'resets'
    return g


for _n, _all in (('reset', False), ('resetAll', True), ('resetToOrigin', True)):
    setattr(ExtendedRoadRunner, _n, _resets(getattr(roadrunner.RoadRunner, _n), _all))


# Structure-changing roadrunner methods drop the cached Rust model (rebuilt on next use)
def _invalidates(f):
    def g(self, *args, **kwargs):
        if f.__name__.startswith('load'):
            self.__dict__.pop('_lz', None)  # a new model replaces the lazy one
        self._load_rr().__dict__.pop('_rs_cache', None)
        return f(self, *args, **kwargs)
    g.__name__, g.__doc__, g._route = f.__name__, f.__doc__, 'invalidates'
    return g


for _n in dir(roadrunner.RoadRunner):
    if _n.startswith(('add', 'remove', 'load', 'clear', 'regenerate', 'setBoundary', 'setConstant',
                      'setHasOnlySubstanceUnits', 'setInitStoichiometryValue', 'setKineticLaw', 'setReversible')):
        setattr(ExtendedRoadRunner, _n, _invalidates(getattr(roadrunner.RoadRunner, _n)))


# Every other roadrunner method needs roadrunner's model: compile it first when the model is lazy.
# These don't (solver settings live in the RoadRunner object and survive the load).
_NO_MODEL = {'getAvailableIntegrators', 'getDiffStepSize', 'setDiffStepSize', 'getExistingIntegratorNames',
             'getExistingSteadyStateSolverNames', 'getFluxThreshold',
             'setFluxThreshold', 'getInstanceCount', 'getInstanceID', 'getIntegrator',
             'getSteadyStateSolver', 'getSteadyStateSolverByName', 'getSteadyStateThreshold',
             'setSteadyStateThreshold', 'integratorExists', 'makeSteadyStateSolver', 'setSteadyStateSolver',
             'steadyStateSolverExists',
             'getCompiler', 'showPlot', 'plotLegend'}


def _loads(name):
    def g(self, *args, **kwargs):
        return getattr(self._super(), name)(*args, **kwargs)
    g.__name__, g.__doc__, g._route = name, getattr(roadrunner.RoadRunner, name).__doc__, 'loads'
    return g


def _public(n, f):
    return inspect.isfunction(f) and (not n.startswith('_') or n == '__getstate__')


for _n, _f in list(vars(roadrunner.RoadRunner).items()):
    if _public(_n, _f) and _n not in vars(ExtendedRoadRunner) and _n not in _NO_MODEL:
        setattr(ExtendedRoadRunner, _n, _loads(_n))


def _routing():
    """How each public roadrunner method is routed: name -> 'override' (defined here: Rust when it can),
    'loads' (compiles a lazy model, then roadrunner), 'invalidates' (drops the Rust model), 'resets'
    (restarts events), 'no-model' (needs no compiled model) or 'unrouted' (roadrunner's own, which
    would run on an uncompiled lazy model: a bug)."""
    out = {}
    for n, f in vars(roadrunner.roadrunner.RoadRunner).items():  # roadrunner.RoadRunner is patched by tellurium
        if _public(n, f):
            g = getattr(ExtendedRoadRunner, n)
            out[n] = 'no-model' if n in _NO_MODEL else getattr(g, '_route', 'unrouted' if g is f else 'override')
    return out


def _lzm(f):  # _LazyModel method: lazily from the Rust state, else roadrunner's model
    def g(self, *args):
        lz = self._r.__dict__.get('_lz')
        if lz is None:
            return getattr(self._r.model, f.__name__)(*args)
        return f(self, lz, *args)
    g.__name__ = f.__name__
    return g


class _LazyModel(object):
    """``r.model`` before roadrunner compiled the model: values from the Rust model's state. Anything
    else compiles the model in roadrunner and is answered by roadrunner's model."""

    def __init__(self, r):
        self._r = r

    def __getattr__(self, name):
        return getattr(self._r._load_rr().model, name)

    def __getitem__(self, key):
        return self._r.getValue(key)

    def __setitem__(self, key, value):
        self._r.setValue(key, value)

    def _vals(self, lz, ids, kind, idx):
        ids = ids if idx is None else [ids[i] for i in idx]
        st, out = lz['st'], self._r._lz_out(lz)
        v = out[lz['o'](ids)]
        if kind is None:
            return v
        V = out[lz['o']([st.comp_of[s] for s in ids])]
        a = np.array([s in st.amount_params for s in ids], dtype=bool)  # value is an amount
        return np.where(a, v / V, v) if kind == 'conc' else np.where(a, v, v * V)

    def _set(self, ids, fmt, args):
        idx, vals = (range(len(ids)), args[0]) if len(args) == 1 else args
        for i, v in zip(idx, vals):
            self._r.setValue(fmt % ids[i], v)

    @_lzm
    def getTime(self, lz):
        return lz['t']

    @_lzm
    def setTime(self, lz, t):
        lz['t'] = float(t)

    @_lzm
    def getNumEvents(self, lz):
        return lz['info']['doc'].getModel().getNumEvents()

    @_lzm
    def getNumRateRules(self, lz):
        return len(lz['info']['rate_rule_ids'])

    @_lzm
    def getNumConservedMoieties(self, lz):  # conservedMoietyAnalysis is off while lazy
        return 0

    @_lzm
    def getConservedMoietyIds(self, lz):
        return []

    @_lzm
    def getConservedMoietyValues(self, lz):
        return np.zeros(0)

    @_lzm
    def getFloatingSpeciesInitConcentrations(self, lz):
        v = self._r._lz_init_values(lz) if lz['reset_ok'] else lz['x0'][lz['tgt']]
        d = dict(zip(lz['tgt_ids'], v))
        c = self.getFloatingSpeciesConcentrations()
        return np.array([d.get(s, x) for s, x in zip(lz['st'].ids, c)])


def _lazy_ids(key, kind=None, fmt=None):  # _LazyModel getters/setters over one id list
    def ids(self, lz):
        return list(getattr(lz['st'], key))

    def num(self, lz):
        return len(getattr(lz['st'], key))

    def get(self, lz, idx=None):
        return self._vals(lz, getattr(lz['st'], key), kind, idx)

    def put(self, lz, *args):
        self._set(getattr(lz['st'], key), fmt, args)
    return ids, num, get, put


for _k, _i, _n, _g, _s, _kind, _fmt in (
        ('ids', 'getFloatingSpeciesIds', 'getNumFloatingSpecies', 'getFloatingSpeciesConcentrations',
         'setFloatingSpeciesConcentrations', 'conc', '[%s]'),
        ('ids', None, None, 'getFloatingSpeciesAmounts', 'setFloatingSpeciesAmounts', 'amt', '%s'),
        ('bids', 'getBoundarySpeciesIds', 'getNumBoundarySpecies', 'getBoundarySpeciesConcentrations',
         'setBoundarySpeciesConcentrations', 'conc', '[%s]'),
        ('bids', None, None, 'getBoundarySpeciesAmounts', 'setBoundarySpeciesAmounts', 'amt', '%s'),
        ('gids', 'getGlobalParameterIds', 'getNumGlobalParameters', 'getGlobalParameterValues',
         'setGlobalParameterValues', None, '%s'),
        ('cids', 'getCompartmentIds', 'getNumCompartments', 'getCompartmentVolumes', 'setCompartmentVolumes', None, '%s'),
        ('rids', 'getReactionIds', 'getNumReactions', 'getReactionRates', None, None, '%s')):
    for _name, _fn in zip((_i, _n, _g, _s), _lazy_ids(_k, _kind, _fmt)):
        if _name:
            _fn.__name__ = _name
            setattr(_LazyModel, _name, _lzm(_fn))
