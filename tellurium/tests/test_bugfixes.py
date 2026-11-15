"""Regression tests for core bug fixes."""
import os
import tempfile
import types

import tellurium as te


def test_null_engine_tiled_and_lowertri_figures():
    old = te.getDefaultPlottingEngine()
    try:
        te.disablePlotting()
        te.newTiledFigure(rows=1, cols=1)
        te.newLowerTriFigure(rows=2, cols=2)
    finally:
        te.setDefaultPlottingEngine(old)


class _Stop(Exception):
    pass


class _FakeSpark(object):
    def addPyFile(self, f):
        pass

    def parallelize(self, samples, n):
        assert len(samples) == n == 3
        raise _Stop()


def test_distributed_sensitivity_analysis_two_value_bounds():
    for allow_log in (False, True):
        model = types.SimpleNamespace(bounds={'k1': [0, 1]}, filename='/x/model.py',
                                      allowLog=allow_log)
        try:
            te.distributed_sensitivity_analysis(_FakeSpark(), model, calculation='avg')
        except _Stop:
            pass
        else:
            raise AssertionError('did not reach Spark')


def test_inline_omex_without_phrasedml(monkeypatch):
    import importlib
    import sys
    import pytest
    from tellurium.teconverters import inline_omex
    monkeypatch.setitem(sys.modules, 'phrasedml', None)  # simulate not installed
    try:
        importlib.reload(inline_omex)
        assert hasattr(inline_omex, 'saveInlineOMEX')
        with pytest.raises(ImportError):
            inline_omex.inlineOmex([])
    finally:
        monkeypatch.undo()
        importlib.reload(inline_omex)


def test_execute_omex_removes_working_dir(monkeypatch):
    from tellurium.teconverters.convert_omex import Omex
    created = []
    real_mkdtemp = tempfile.mkdtemp

    def mkdtemp(*a, **kw):
        d = real_mkdtemp(*a, **kw)
        created.append(d)
        return d

    monkeypatch.setattr(tempfile, 'mkdtemp', mkdtemp)
    monkeypatch.setattr(te, 'executeSEDML', lambda *a, **kw: None)
    Omex(description=None, creator=None).executeOmex()
    assert created and not any(os.path.exists(d) for d in created)
