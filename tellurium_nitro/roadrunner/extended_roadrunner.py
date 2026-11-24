"""
tellurium's RoadRunner: libroadrunner's Python API on the Rust engine (roadrunner_nitro.RoadRunner).
"""

from __future__ import print_function, division, absolute_import
import copy
import os
import re
import tempfile
import warnings
from xml.sax.saxutils import quoteattr

import numpy as np
import roadrunner_nitro

INTEGRATORS = ['cvode', 'gillespie', 'rk4', 'rk45', 'euler']
STEADY_STATE_SOLVERS = ['nleq1', 'nleq2', 'newton', 'newton_linesearch']


def turbo(enable=None):
    """Kept for compatibility. The Rust engine is tellurium_nitro's only engine (libroadrunner is no longer
    a fallback), so turbo is always on: True and 'strict' are accepted, False warns and changes nothing.

    :returns: True
    """
    if enable is not None and not enable:
        warnings.warn("te.turbo(False): there is no other engine than the Rust engine; it stays on.", stacklevel=2)
    return True


# ---------------------------------------------------------------------
# NamedArray, Logger, Config, SelectionRecord (libroadrunner's helper types)
# ---------------------------------------------------------------------
class NamedArray(np.ndarray):
    """numpy array with ``rownames`` / ``colnames``, printed as libroadrunner prints its NamedArray.
    ``a['S1']``: the column (or row) of that name."""

    def __new__(cls, a, rownames=(), colnames=()):
        obj = np.asarray(a).view(cls)
        obj.rownames, obj.colnames = list(rownames), list(colnames)
        return obj

    def __array_finalize__(self, obj):
        self.rownames = list(getattr(obj, 'rownames', ()))
        self.colnames = list(getattr(obj, 'colnames', ()))

    def __reduce__(self):
        f, args, state = super(NamedArray, self).__reduce__()
        return f, args, (state, self.rownames, self.colnames)

    def __setstate__(self, state):
        super(NamedArray, self).__setstate__(state[0])
        self.rownames, self.colnames = state[1], state[2]

    def __getitem__(self, key):
        if isinstance(key, str):
            a = np.asarray(self)
            if key in self.colnames:
                return a[..., self.colnames.index(key)]
            if key in self.rownames:
                return a[self.rownames.index(key)]
            raise KeyError(key)
        return super(NamedArray, self).__getitem__(key)

    def __str__(self):
        a = np.asarray(self)
        if a.ndim != 2 or not len(a) or len(self.colnames) != a.shape[1] or \
                (self.rownames and len(self.rownames) != a.shape[0]):
            return str(a)
        cells = [[format(v, 'g') for v in row] for row in a]
        w = [max([len(n)] + [len(r[j]) for r in cells]) + 1 for j, n in enumerate(self.colnames)]
        rn = max([len(n) for n in self.rownames] or [0])
        out = [' ' * rn + '   ' + ','.join(n.rjust(x) for n, x in zip(self.colnames, w))]
        for i, r in enumerate(cells):
            row = ','.join(c.rjust(x) for c, x in zip(r, w))
            lead = self.rownames[i].ljust(rn) if rn else ''
            out.append(lead + (' [[' if i == 0 else '  [') + row + (']]' if i == len(cells) - 1 else '],'))
        return '\n'.join(out) + '\n'

    __repr__ = __str__


def _named(t):
    """(array, row ids, column ids) from roadrunner_nitro -> NamedArray."""
    return NamedArray(t[0], t[1], t[2])


class Logger(object):
    """libroadrunner's Logger API (levels as libroadrunner numbers them), on Python's logging."""
    LOG_CURRENT, LOG_FATAL, LOG_CRITICAL, LOG_ERROR, LOG_WARNING, LOG_NOTICE, LOG_INFORMATION, LOG_DEBUG, LOG_TRACE = range(9)
    _level = LOG_NOTICE

    @classmethod
    def setLevel(cls, level=LOG_CURRENT):
        if level != cls.LOG_CURRENT:
            cls._level = level

    @classmethod
    def getLevel(cls):
        return cls._level

    @classmethod
    def log(cls, level, msg):
        if level <= cls._level:
            import logging
            logging.getLogger('tellurium_nitro').log(max(10, 60 - 10 * level), msg)


class Config(object):
    """libroadrunner's Config: setValue/getValue on a key; the Rust engine reads none of them."""
    _values = {}

    @classmethod
    def setValue(cls, key, value):
        cls._values[key] = value

    @classmethod
    def getValue(cls, key):
        return cls._values.get(key)


class SelectionRecord(object):
    """libroadrunner's selection type flags (the ones reset(options) takes)."""
    TIME, CONCENTRATION, AMOUNT, RATE, BOUNDARY, FLOATING, _COMPARTMENT, _GLOBAL_PARAMETER = 1, 2, 4, 8, 16, 32, 64, 128
    REACTION, INITIAL, CURRENT, UNSCALED, ELASTICITY, CONTROL = 256, 512, 1024, 2048, 4096, 8192
    ALL = -1


# ---------------------------------------------------------------------
# Solver settings objects: r.integrator, r.steadyStateSolver
# ---------------------------------------------------------------------
class Solver(object):
    """An integrator or steady-state solver: settings as attributes, getValue/setValue/getSettings.
    Without ``name`` it is the current one; with ``name``, that solver (current or not)."""

    def __init__(self, rr, ss, name=None):
        object.__setattr__(self, '_rr', rr)
        object.__setattr__(self, '_ss', ss)
        object.__setattr__(self, '_name', name)

    def _call(self, f, *args):
        rr, ss = self._rr, self._ss
        cur = rr.getSteadyStateSolverName() if ss else rr.getIntegratorName()
        use = rr.setSteadyStateSolver if ss else rr.setIntegrator
        if self._name is None or self._name == cur:
            return f(*args)
        use(self._name)  # settings of the other solvers are kept: switch, act, switch back
        try:
            return f(*args)
        finally:
            use(cur)

    def getName(self):
        if self._name is not None:
            return self._name
        return self._rr.getSteadyStateSolverName() if self._ss else self._rr.getIntegratorName()

    def getSettings(self):
        return list(self._call(self._rr.getSteadyStateSolverSettings if self._ss else self._rr.getIntegratorSettings))

    def getValue(self, key):
        return self._call(self._rr.getSteadyStateSolverSetting if self._ss else self._rr.getIntegratorSetting, key)

    def setValue(self, key, value):
        if isinstance(value, np.generic):
            value = value.item()
        self._call(self._rr.setSteadyStateSolverSetting if self._ss else self._rr.setIntegratorSetting, key, value)

    def getSettingsMap(self):
        return {k: self.getValue(k) for k in self.getSettings()}

    def getSetting(self, key):
        return self.getValue(key)

    def setSetting(self, key, value):
        self.setValue(key, value)

    def __getattr__(self, key):
        if key.startswith('_'):
            raise AttributeError(key)
        try:
            return self.getValue(key)
        except ValueError:
            raise AttributeError("%s has no setting '%s'" % (self.getName(), key))

    def __setattr__(self, key, value):
        self.setValue(key, value)

    def __dir__(self):
        return sorted(set(dir(type(self))) | set(self.getSettings()))

    def __repr__(self):
        s = self.getSettingsMap()
        w = max([len(k) for k in s] or [0])
        kind = 'SteadyStateSolver' if self._ss else 'Integrator'
        return '< roadrunner.%s() >\n  name: %s\n  settings:\n%s' % (kind, self.getName(), ''.join(
            '  %s: %s\n' % (k.rjust(w + 4), str(v).lower() if isinstance(v, bool) else '%g' % v if isinstance(v, float) else v) for k, v in s.items()))

    __str__ = __repr__


class SensitivitySolver(Solver):
    """r.getSensitivitySolver(): the tolerances timeSeriesSensitivities integrates with ('forward': exact
    symbolic derivatives, one stiff integration per parameter)."""

    def __init__(self):
        object.__setattr__(self, '_s', {'relative_tolerance': 1e-6, 'absolute_tolerance': 1e-12})

    def getName(self):
        return 'forward'

    def getSettings(self):
        return list(self._s)

    def getValue(self, key):
        if key not in self._s:
            raise ValueError("forward: no setting '%s'" % key)
        return self._s[key]

    def setValue(self, key, value):
        self.getValue(key)
        self._s[key] = float(value)


# ---------------------------------------------------------------------
# r.model: libroadrunner's ExecutableModel API
# ---------------------------------------------------------------------
def _sbml_ids(sbml, tag, attr):
    """Ids from the SBML text: <tag ... attr="x"> (rule variables, initial assignment symbols, event ids)."""
    return re.findall(r'<(?:\w+:)?%s\b[^>]*?\b%s\s*=\s*"([^"]*)"' % (tag, attr), sbml)


class ExecutableModel(object):
    """``r.model``: ids and values of the loaded model, as libroadrunner's ExecutableModel."""

    def __init__(self, r):
        self._r = r

    def _get(self, ids, fmt, idx=None):
        ids = ids if idx is None else [ids[i] for i in np.atleast_1d(idx)]
        return np.array(self._r._rr.getValues([fmt % i for i in ids]), dtype=float)

    def _set(self, ids, fmt, args):
        idx, vals = (range(len(ids)), args[0]) if len(args) == 1 else args
        for i, v in zip(np.atleast_1d(idx), np.atleast_1d(vals)):
            self._r._rr.setValue(fmt % ids[i], float(v))

    def __getitem__(self, key):
        return self._r.getValue(key)

    def __setitem__(self, key, value):
        self._r.setValue(key, value)

    def keys(self):
        return [s for _, l in self._r._rr.getAvailableTimeCourseSymbols() for s in l]

    def values(self):
        return [self[k] for k in self.keys()]

    def items(self):
        return list(zip(self.keys(), self.values()))

    def getTime(self):
        return self._r._rr.getCurrentTime()

    def setTime(self, t):
        self._r._rr.setValue('time', float(t))

    def reset(self):
        self._r.reset()

    def getModelName(self):
        return self._r.getModelName()

    def getInfo(self):
        return self._r.getInfo()

    def getNumEvents(self):
        return len(self.getEventIds())

    def getEventIds(self):
        return _sbml_ids(self._r._rr.getSBML(), 'event', 'id')

    def getAssignmentRuleIds(self):
        return _sbml_ids(self._r._rr.getSBML(), 'assignmentRule', 'variable')

    def getInitialAssignmentIds(self):
        return _sbml_ids(self._r._rr.getSBML(), 'initialAssignment', 'symbol')

    def getRateRuleIds(self):
        return list(self._r._rr.getRateRuleIds())

    def getNumRateRules(self):
        return len(self.getRateRuleIds())

    def getRateRuleRates(self):
        return self._get(self.getRateRuleIds(), "%s'")

    def getConservedMoietyIds(self):
        return list(self._r._rr.getConservedMoietyIds())

    def getNumConservedMoieties(self):
        return len(self.getConservedMoietyIds())

    def getConservedMoietyValues(self, idx=None):
        return self._get(self.getConservedMoietyIds(), '%s', idx)

    def setConservedMoietyValues(self, *args):
        self._set(self.getConservedMoietyIds(), '%s', args)

    def getFloatingSpeciesConcentrationIds(self):
        return ['[%s]' % s for s in self._r._rr.getFloatingSpeciesIds()]

    def getBoundarySpeciesConcentrationIds(self):
        return ['[%s]' % s for s in self._r._rr.getBoundarySpeciesIds()]

    def getFloatingSpeciesInitConcentrationIds(self):
        return ['init([%s])' % s for s in self._r._rr.getFloatingSpeciesIds()]

    def getFloatingSpeciesInitAmountIds(self):
        return ['init(%s)' % s for s in self._r._rr.getFloatingSpeciesIds()]

    def getConservedMoietyId(self, i):
        return self.getConservedMoietyIds()[i]

    def getEventId(self, i):
        return self.getEventIds()[i]

    def getFloatingSpeciesAmountRate(self, i, reactionRates=None):
        """Amount rate of floating species i from the reactions: N[i] . reactionRates (default: the current rates)."""
        n, rows, cols = self._r._rr.getFullStoichiometryMatrix()
        v = self.getReactionRates() if reactionRates is None else np.asarray(reactionRates)
        ids = self._r._rr.getReactionIds()
        return float(np.asarray(n)[rows.index(self._r._rr.getFloatingSpeciesIds()[i]), [cols.index(r) for r in ids]] @ v)

    def _compartment_index(self, sid):
        m = re.search(r'<(?:\w+:)?species\b[^>]*?\bid\s*=\s*"%s"[^>]*' % re.escape(sid), self._r._rr.getSBML())
        c = re.search(r'\bcompartment\s*=\s*"([^"]*)"', m.group(0)).group(1)
        return self._r._rr.getCompartmentIds().index(c)

    def getCompartmentIndexForFloatingSpecies(self, i):
        return self._compartment_index(self._r._rr.getFloatingSpeciesIds()[i])

    def getCompartmentIndexForBoundarySpecies(self, i):
        return self._compartment_index(self._r._rr.getBoundarySpeciesIds()[i])

    def getStoichiometry(self, speciesIndex, reactionIndex):
        n, rows, cols = self._r._rr.getFullStoichiometryMatrix()
        rr = self._r._rr
        return float(np.asarray(n)[rows.index(rr.getFloatingSpeciesIds()[speciesIndex]), cols.index(rr.getReactionIds()[reactionIndex])])

    def getStateVectorId(self, i):
        return list(self._r._rr.getRateRuleIds() + self._r._rr.getFloatingSpeciesIds())[i]


def _model_group(ids, getters, setters):
    """ExecutableModel ids / count / value getters (name -> selection format) / setters for one id list."""
    def ids_(self):
        return list(getattr(self._r._rr, ids)())

    def num(self):
        return len(getattr(self._r._rr, ids)())

    def getter(fmt):
        return lambda self, idx=None: self._get(getattr(self._r._rr, ids)(), fmt, idx)

    def setter(fmt):
        return lambda self, *args: self._set(getattr(self._r._rr, ids)(), fmt, args)
    group = ids[3:-3]  # FloatingSpecies, GlobalParameter, ...
    out = {ids: ids_, 'getNum' + group + ('' if group.endswith('Species') else 's'): num,
           'get%sId' % group: lambda self, i: getattr(self._r._rr, ids)()[i]}
    out.update({n: getter(f) for n, f in getters.items()})
    out.update({n: setter(f) for n, f in setters.items()})
    return out


for _ids, _get, _put in (
        ('getFloatingSpeciesIds', {'getFloatingSpeciesConcentrations': '[%s]', 'getFloatingSpeciesAmounts': '%s',
                                   'getFloatingSpeciesInitConcentrations': 'init([%s])',
                                   'getFloatingSpeciesInitAmounts': 'init(%s)',
                                   'getFloatingSpeciesConcentrationRates': "[%s]'",
                                   'getFloatingSpeciesAmountRates': "%s'"},
         {'setFloatingSpeciesConcentrations': '[%s]', 'setFloatingSpeciesAmounts': '%s',
          'setFloatingSpeciesInitConcentrations': 'init([%s])', 'setFloatingSpeciesInitAmounts': 'init(%s)'}),
        ('getBoundarySpeciesIds', {'getBoundarySpeciesConcentrations': '[%s]', 'getBoundarySpeciesAmounts': '%s',
                                   'getBoundarySpeciesInitConcentrations': 'init([%s])',
                                   'getBoundarySpeciesInitAmounts': 'init(%s)'},
         {'setBoundarySpeciesConcentrations': '[%s]', 'setBoundarySpeciesAmounts': '%s',
          'setBoundarySpeciesInitConcentrations': 'init([%s])', 'setBoundarySpeciesInitAmounts': 'init(%s)'}),
        ('getGlobalParameterIds', {'getGlobalParameterValues': '%s', 'getGlobalParameterInitValues': 'init(%s)'},
         {'setGlobalParameterValues': '%s', 'setGlobalParameterInitValues': 'init(%s)'}),
        ('getCompartmentIds', {'getCompartmentVolumes': '%s', 'getCompartmentInitVolumes': 'init(%s)'},
         {'setCompartmentVolumes': '%s', 'setCompartmentInitVolumes': 'init(%s)'}),
        ('getReactionIds', {'getReactionRates': '%s'}, {})):
    for _n, _f in _model_group(_ids, _get, _put).items():
        _f.__name__ = _n
        setattr(ExecutableModel, _n, _f)
ExecutableModel.getStateVectorRate = ExecutableModel.getFloatingSpeciesAmountRates
ExecutableModel.getRateRuleSymbols = ExecutableModel.getRateRuleIds
ExecutableModel.getAllTimeCourseComponentIds = ExecutableModel.keys
ExecutableModel.getValue = ExecutableModel.__getitem__
ExecutableModel.setValue = ExecutableModel.__setitem__
ExecutableModel.iterkeys, ExecutableModel.itervalues, ExecutableModel.iteritems = \
    ExecutableModel.keys, ExecutableModel.values, ExecutableModel.items


# ---------------------------------------------------------------------
# Model editing: libroadrunner's method -> (roadrunner_nitro.sbml_edit op, argument defaults after the ids)
# ---------------------------------------------------------------------
_EDITS = {
    'addSpeciesAmount': ('add_species_amount', 2, (0.0, False, False, '')),
    'addSpeciesConcentration': ('add_species_concentration', 2, (0.0, False, False, '')),
    'removeSpecies': ('remove_species', 1, ()),
    'setBoundary': ('set_boundary', 2, ()),
    'setHasOnlySubstanceUnits': ('set_has_only_substance_units', 2, ()),
    'setConstant': ('set_constant', 2, ()),
    'removeReaction': ('remove_reaction', 1, (False,)),
    'setReversible': ('set_reversible', 2, ()),
    'setKineticLaw': ('set_kinetic_law', 2, ()),
    'addParameter': ('add_parameter', 2, ()),
    'removeParameter': ('remove_parameter', 1, ()),
    'addCompartment': ('add_compartment', 2, ()),
    'removeCompartment': ('remove_compartment', 1, ()),
    'addInitialAssignment': ('add_initial_assignment', 2, ()),
    'addAssignmentRule': ('add_assignment_rule', 2, ()),
    'addRateRule': ('add_rate_rule', 2, ()),
    'addEvent': ('add_event', 3, ()),
    'addTrigger': ('add_trigger', 2, ()),
    'setPersistent': ('set_persistent', 2, ()),
    'setTriggerInitialValue': ('set_trigger_initial_value', 2, ()),
    'addPriority': ('add_priority', 2, ()),
    'addDelay': ('add_delay', 2, ()),
    'addEventAssignment': ('add_event_assignment', 3, ()),
    'removeEventAssignments': ('remove_event_assignments', 2, ()),
    'removeEvent': ('remove_event', 1, ()),
}


def _edit_method(name, op, n, defaults):
    def f(self, *args, **kwargs):
        regen = kwargs.pop('forceRegenerate', True)
        args = list(args) + [None] * (n + len(defaults) - len(args))
        if len(args) > n + len(defaults):  # forceRegenerate given positionally
            regen = args.pop()
        for i, d in enumerate(defaults):
            if args[n + i] is None:
                args[n + i] = d
        names = f._kw
        for k, v in kwargs.items():
            args[names.index(k)] = v
        self._edit(op, args, regen)
    f.__name__ = name
    f._kw = {'addSpeciesAmount': ['sid', 'compartment', 'initAmount', 'hasOnlySubstanceUnits', 'boundaryCondition', 'substanceUnits'],
             'addSpeciesConcentration': ['sid', 'compartment', 'initConcentration', 'hasOnlySubstanceUnits',
                                         'boundaryCondition', 'substanceUnits'],
             'removeReaction': ['rid', 'deleteUnusedParameters']}.get(name, [])
    f.__doc__ = "libroadrunner's %s (SBML edited in Rust; forceRegenerate=True recompiles keeping the current values)." % name
    return f


def _set_attr(sbml, tag, id_, attr, value, id_attr='id'):
    """SBML text with attribute ``attr`` of the <tag id_attr="id_"> start tag (id_ None: the first <tag>) set to ``value``."""
    sel = r'[^>]*?' if id_ is None else r'[^>]*?\b%s\s*=\s*"%s"[^>]*?' % (id_attr, re.escape(id_))
    pat = re.compile(r'<((?:\w+:)?%s)\b(%s)(/?)>' % (tag, sel))
    m = pat.search(sbml)
    if m is None:
        raise ValueError("no %s '%s' in the model" % (tag, id_))
    body = re.sub(r'\s%s\s*=\s*"[^"]*"' % attr, '', m.group(2)) + ' %s=%s' % (attr, quoteattr(str(value)))
    return sbml[:m.start()] + '<%s%s%s>' % (m.group(1), body, m.group(3)) + sbml[m.end():]


# ---------------------------------------------------------------------
# RoadRunner
# ---------------------------------------------------------------------
class ExtendedRoadRunner(object):
    """libroadrunner's RoadRunner API on the Rust engine, plus tellurium's helpers (plot, export, ...).

    RoadRunner(sbml=None): SBML text or a file path."""

    def __init__(self, sbml=None, options=None):
        # options (libroadrunner's LoadSBMLOptions): the Rust engine has no load options
        object.__setattr__(self, '_rr', roadrunner_nitro.RoadRunner())
        object.__setattr__(self, '_data', None)
        object.__setattr__(self, '_last_run', None)
        object.__setattr__(self, '_attrs', None)
        object.__setattr__(self, 'model', ExecutableModel(self))
        object.__setattr__(self, '_sens', SensitivitySolver())
        if sbml is not None:
            self.load(sbml)

    # ------------------------------------------------------------ loading, SBML
    def load(self, sbml):
        """Load an SBML model from text or a file path."""
        if isinstance(sbml, bytes):
            sbml = sbml.decode('utf-8')
        if sbml.startswith('file://'):
            sbml = sbml[len('file://'):]
        self._rr.load(sbml)
        self._changed()

    def _changed(self):  # the model's ids may have changed
        object.__setattr__(self, '_attrs', None)

    def clearModel(self):
        self._rr.clearModel()
        self._changed()

    def isModelLoaded(self):
        return self._rr.isModelLoaded()

    def getModel(self):
        return self.model

    def getSBML(self, level=0, version=0):
        """The SBML as loaded (with model edits)."""
        if level or version:
            raise ValueError("getSBML: converting to another SBML level/version is not supported")
        return self._rr.getSBML()

    def getCurrentSBML(self, level=0, version=0):
        """The SBML with the current values as initial values."""
        if level or version:
            raise ValueError("getCurrentSBML: converting to another SBML level/version is not supported")
        return self._rr.getCurrentSBML()

    def getParamPromotedSBML(self, sbml=None):
        """SBML with local parameters promoted to global ones (of ``sbml``, default the current model)."""
        return roadrunner_nitro.sbml_edit('get_param_promoted_sbml', self.getSBML() if sbml is None else sbml)

    def validateCurrentSBML(self):
        return roadrunner_nitro.sbml_edit('validate_current_sbml', self.getSBML())

    def getModelName(self):
        m = re.search(r'<(?:\w+:)?model\b[^>]*?\bname\s*=\s*"([^"]*)"', self._rr.getSBML())
        return m.group(1) if m else self.getModelId()

    def getModelId(self):
        m = re.search(r'<(?:\w+:)?model\b[^>]*?\bid\s*=\s*"([^"]*)"', self._rr.getSBML())
        return m.group(1) if m else ''

    def getInfo(self):
        rr = self._rr
        return '\n'.join(['<roadrunner.RoadRunner() >', 'engine: Rust (roadrunner_nitro)',
                          'model loaded: %s' % rr.isModelLoaded(), 'modelName: %s' % self.getModelName(),
                          'integrator: %s' % rr.getIntegratorName(), 'conservedMoietyAnalysis: %s' % rr.conservedMoietyAnalysis,
                          'floating species: %d' % len(rr.getFloatingSpeciesIds()),
                          'reactions: %d' % len(rr.getReactionIds())])

    def getExtendedVersionInfo(self):
        return 'roadrunner_nitro (Rust engine)'

    # ------------------------------------------------------------ values
    def getValue(self, sel):
        return self._rr.getValue(sel)

    def setValue(self, sel, value):
        self._rr.setValue(sel, float(value))

    def setValues(self, keysOrDict, values=None):
        items = keysOrDict.items() if values is None else zip(keysOrDict, values)
        for k, v in items:
            self.setValue(k, v)

    def __getitem__(self, sel):
        if isinstance(sel, (list, tuple)):
            return np.array(self._rr.getValues(list(sel)))
        return self._rr.getValue(sel)

    def __setitem__(self, sel, value):
        self.setValue(sel, value)

    def keys(self):
        return self.model.keys()

    def values(self):
        return self.model.values()

    def items(self):
        return self.model.items()

    iterkeys, itervalues, iteritems = keys, values, items

    def resetSelectionLists(self):
        """Selections back to the defaults: time, floating species concentrations, rate-rule values."""
        fs = self._rr.getFloatingSpeciesIds()
        self.timeCourseSelections = ['time'] + ['[%s]' % s for s in fs] + [v for v in self._rr.getRateRuleIds() if v not in fs]
        self.steadyStateSelections = ['[%s]' % s for s in fs]

    def getHasOnlySubstanceUnits(self, sid):
        m = re.search(r'<(?:\w+:)?species\b[^>]*?\bid\s*=\s*"%s"[^>]*' % re.escape(sid), self._rr.getSBML())
        if m is None:
            raise ValueError("no species '%s'" % sid)
        return re.search(r'\bhasOnlySubstanceUnits\s*=\s*"(true|1)"', m.group(0)) is not None

    def getSelectedValues(self):
        return np.array(self._rr.getValues(self.timeCourseSelections))

    def _attr_map(self):
        """Model ids as attributes: r.S1 (concentration), r.S1_amt, r.S1_conc, parameters, compartments, rates."""
        a = self._attrs
        if a is None:
            rr, a = self._rr, {}
            for s in rr.getFloatingSpeciesIds():
                a.update({s: '[%s]' % s, s + '_conc': '[%s]' % s, s + '_amt': s})
            a.update({s: '[%s]' % s for s in rr.getBoundarySpeciesIds()})
            a.update({k: k for k in rr.getGlobalParameterIds() + rr.getCompartmentIds() + rr.getReactionIds()})
            object.__setattr__(self, '_attrs', a)
        return a

    def __getattr__(self, name):
        if not name.startswith('_'):
            sel = self._attr_map().get(name)
            if sel is not None:
                return self._rr.getValue(sel)
        raise AttributeError("'RoadRunner' object has no attribute '%s'" % name)

    def __setattr__(self, name, value):
        if not name.startswith('_') and not hasattr(type(self), name):
            sel = self._attr_map().get(name)
            if sel is not None:
                return self.setValue(sel, value)
        object.__setattr__(self, name, value)

    def __dir__(self):
        return sorted(set(dir(type(self))) | set(self._attr_map()))

    def getCurrentTime(self):
        return self._rr.getCurrentTime()

    def _ran(self, method):
        object.__setattr__(self, '_last_run', dict(method=method, engine='rust', reason=None))

    @property
    def lastRun(self):
        """The last simulate/steadyState/MCA-type call: dict(method, engine='rust', reason=None);
        None before any such call. (Kept from the two-engine era: the engine is always Rust.)"""
        r = self._last_run
        return None if r is None else dict(r)

    # ------------------------------------------------------------ selections
    timeCourseSelections = selections = property(
        lambda self: list(self._rr.timeCourseSelections),
        lambda self, s: setattr(self._rr, 'timeCourseSelections', list(s)))
    steadyStateSelections = property(
        lambda self: list(self._rr.steadyStateSelections),
        lambda self, s: setattr(self._rr, 'steadyStateSelections', list(s)))

    def getSteadyStateSelectionStrings(self):
        return self.steadyStateSelections

    @property
    def conservedMoietyAnalysis(self):
        return self._rr.conservedMoietyAnalysis

    @conservedMoietyAnalysis.setter
    def conservedMoietyAnalysis(self, on):
        self._rr.conservedMoietyAnalysis = bool(on)
        self._changed()

    # ------------------------------------------------------------ solvers
    def getIntegrator(self):
        return Solver(self._rr, False)

    def setIntegrator(self, name):
        self._rr.setIntegrator(name if isinstance(name, str) else name.getName())

    integrator = property(getIntegrator, setIntegrator)

    def getIntegratorByName(self, name):
        if name not in INTEGRATORS:
            raise ValueError("no integrator '%s' (one of %s)" % (name, INTEGRATORS))
        return Solver(self._rr, False, name)

    makeIntegrator = getIntegratorByName

    def setIntegratorSetting(self, integratorName, settingName, value):
        self.getIntegratorByName(integratorName).setValue(settingName, value)

    def integratorExists(self, name):
        return name in INTEGRATORS

    def getExistingIntegratorNames(self):
        return list(INTEGRATORS)

    getAvailableIntegrators = getRegisteredIntegratorNames = getExistingIntegratorNames

    def getSteadyStateSolver(self):
        return Solver(self._rr, True)

    def setSteadyStateSolver(self, name):
        self._rr.setSteadyStateSolver(name if isinstance(name, str) else name.getName())

    steadyStateSolver = property(getSteadyStateSolver, setSteadyStateSolver)

    def getSteadyStateSolverByName(self, name):
        if name not in STEADY_STATE_SOLVERS:
            raise ValueError("no steady-state solver '%s' (one of %s)" % (name, STEADY_STATE_SOLVERS))
        return Solver(self._rr, True, name)

    makeSteadyStateSolver = getSteadyStateSolverByName

    def steadyStateSolverExists(self, name):
        return name in STEADY_STATE_SOLVERS

    def getExistingSteadyStateSolverNames(self):
        return list(STEADY_STATE_SOLVERS)

    getRegisteredSteadyStateSolverNames = getExistingSteadyStateSolverNames

    def getSensitivitySolver(self):
        return self._sens

    def setSensitivitySolver(self, name):
        self.getSensitivitySolverByName(name if isinstance(name, str) else name.getName())

    def getSensitivitySolverByName(self, name):
        if name != 'forward':
            raise ValueError("no sensitivity solver '%s' (one of ['forward'])" % name)
        return self._sens

    makeSensitivitySolver = getSensitivitySolverByName

    def sensitivitySolverExists(self, name):
        return name == 'forward'

    def getExistingSensitivitySolverNames(self):
        return ['forward']

    getRegisteredSensitivitySolverNames = getExistingSensitivitySolverNames

    def getDiffStepSize(self):
        return self._rr.diffStepSize

    def setDiffStepSize(self, value):
        self._rr.diffStepSize = float(value)

    diffstep = property(getDiffStepSize, setDiffStepSize)

    def getFluxThreshold(self):
        return self._rr.getFluxThreshold()

    def setFluxThreshold(self, value):
        self._rr.setFluxThreshold(float(value))

    def getSeed(self):
        return self._rr.getSeed()

    def setSeed(self, seed, resetModel=True):
        self._rr.setSeed(int(seed))
        if resetModel:
            self.reset()

    # ------------------------------------------------------------ time courses
    def simulate(self, start=None, end=None, points=None, selections=None, output_file=None, steps=None, times=None):
        """Time course from ``start`` to ``end`` at ``points`` points (or ``steps`` + 1, or at ``times``)
        with the current integrator; returns a NamedArray (columns: the selections, which become
        timeCourseSelections). The model is left at the end. ``output_file``: also written as CSV."""
        if isinstance(points, (list, tuple)):  # simulate(start, end, selections)
            points, selections = None, points
        if steps is not None:
            points = int(steps) + 1
        rr = self._rr
        sel = None if selections is None else list(selections)
        # Explicit output points request a grid even when the integrator uses variable steps.
        vss = rr.getIntegratorName() in ('cvode', 'gillespie', 'rk45') and (points is not None or times is not None) and \
            rr.getIntegratorSetting('variable_step_size')
        if vss:
            rr.setIntegratorSetting('variable_step_size', False)
        try:
            if times is not None:
                y = rr.simulateTimes([float(t) for t in times], sel)
            else:
                # as libroadrunner: start 0, end 5, 51 points unless given (each call)
                y = rr.simulate(0.0 if start is None else float(start), 5.0 if end is None else float(end),
                                51 if points is None else int(points), sel)
        finally:
            if vss:
                rr.setIntegratorSetting('variable_step_size', True)
        out = NamedArray(y, colnames=rr.timeCourseSelections)
        object.__setattr__(self, '_data', out)
        self._ran('simulate')
        if output_file:
            np.savetxt(output_file, out, delimiter=',', header=','.join(out.colnames), comments='')
        return out

    def getSimulationData(self):
        return self._data

    def oneStep(self, currentTime, stepSize, reset=True):
        self._ran('oneStep')
        return self._rr.oneStep(float(currentTime), float(stepSize))

    def simulateBatch(self, params, start, end, points, selections=None, threads=None):
        """Many time courses at once on all CPU cores, one per parameter set.

        :param params: dict ``{'k1': [...], 'k2': [...]}`` of equal-length sequences (zipped, not a
            grid), or a pandas DataFrame. Ids can be parameters, compartments, boundary or floating
            species (concentrations), as in ``r[id] = value``.
        :param start, end, points: as in ``simulate``.
        :param selections: output columns (default: timeCourseSelections).
        :param threads: CPU threads (default: all cores).
        :returns: array (n sets, points, len(selections)). Every run starts from the current state;
            the model itself is not changed.
        """
        from ..analysis.rrs_bridge import from_sbml
        if hasattr(params, 'to_dict'):
            params = params.to_dict('list')
        names = list(params)
        vals = np.array([np.atleast_1d(np.asarray(params[k], dtype=float)) for k in names]).T.reshape(-1, len(names))
        sel = list(self.timeCourseSelections if selections is None else selections)
        times = np.linspace(start, end, points)
        sbml = self.getCurrentSBML()
        info = {}
        m = from_sbml(sbml, info)
        comp_of, amt = info['comp_of'], set(info['amount_params'])
        outs, cols = [], []
        for s in sel:  # selection -> (output column, column of the volume to multiply / divide by)
            br = s[:1] == '[' and s[-1:] == ']'
            v = s[1:-1] if br else s
            if v in comp_of:
                f = 0 if br != (v in amt) else (1 if not br else -1)  # amount wanted from conc: x V
                cols.append((v, comp_of[v] if f else None, f))
                outs += [v] + ([comp_of[v]] if f else [])
            else:
                cols.append((s, None, 0))
                outs.append(s)
        outs = list(dict.fromkeys(outs))
        ok = set(outs) <= set(m.output_ids) and all(k.strip('[]') in set(m.state_ids) | set(m.output_ids) for k in names)
        out = np.full((len(vals), points, len(sel)), np.nan)
        if ok and len(vals):
            ig = self.integrator
            ic = dict(zip(self.getFloatingSpeciesIds(), self.model.getFloatingSpeciesInitConcentrations()))
            asc = [abs(ic[s]) if abs(ic.get(s, 0.0)) > 1e-30 else 1.0 for s in m.state_ids]
            Y = np.asarray(m.simulate_batch([k.strip('[]') for k in names], vals, times, threads=threads,
                                            rtol=float(ig.relative_tolerance), atol=float(ig.absolute_tolerance),
                                            outputs=outs, atol_scale=asc))
            o = {n: i for i, n in enumerate(outs)}
            for j, (c, vol, f) in enumerate(cols):
                y = Y[:, :, o[c]]
                out[:, :, j] = y * Y[:, :, o[vol]] if f == 1 else y / Y[:, :, o[vol]] if f == -1 else y
        # selections the batch engine has no column for, and failed runs: one simulate per set
        for i in (range(len(vals)) if not ok else np.flatnonzero(~np.isfinite(out).all(axis=(1, 2)))):
            out[i] = self._one_run(sbml, names, vals[i], start, end, points, sel)
        self._ran('simulateBatch')
        return out

    def _one_run(self, sbml, names, v, start, end, points, sel):
        """One simulateBatch run on a copy of the current model (this one unchanged)."""
        r = ExtendedRoadRunner(sbml)
        r.setIntegrator(self.integrator.getName())
        for k, x in self.integrator.getSettingsMap().items():
            r.integrator.setValue(k, x)
        sp = set(r.getFloatingSpeciesIds()) | set(r.getBoundarySpeciesIds())
        for k, x in zip(names, v):
            r.setValue('[%s]' % k if k in sp else k, float(x))  # species ids are concentrations
        return np.asarray(r.simulate(start, end, points, sel))

    def timeSeriesSensitivities(self, start, stop, num, params=(), species=(), k=0):
        """Forward sensitivities d species / d parameter along a time course (exact symbolic derivatives).

        :returns: (times, sensitivities[time, parameter, species] in amounts, parameter ids, species ids);
            params default: the global parameters; species: rate-rule states, then the (independent)
            floating species. The model is left at ``stop``."""
        from ..analysis.rrs_bridge import from_sbml
        if k not in (0, 1):
            raise ValueError('timeSeriesSensitivities: k must be 0 or 1')
        info = {}
        m = from_sbml(self.getCurrentSBML(), info)
        states, fids = list(m.state_ids), self.getFloatingSpeciesIds()
        cma = self.conservedMoietyAnalysis
        cs = self.getConservedMoietyIds() if cma else []
        rr_ids = [v for v in self.getRateRuleIds() if v in states and v not in fids]
        base = rr_ids + (self.getIndependentFloatingSpeciesIds() if cma else fids)
        params = list(params) or self.getGlobalParameterIds()
        species = [s for s in base if not species or s in species]
        zero = set(rr_ids)
        rp = [p for p in params if p not in cs and p not in zero]
        dep = self.getDependentFloatingSpeciesIds() if any(p in cs for p in params) else []
        times = np.linspace(start, stop, int(num))
        X, S = roadrunner_nitro.time_sensitivities(m, [], [], times, rp + dep, order=k, rtol=self._sens.relative_tolerance,
                                                   atol=self._sens.absolute_tolerance)
        vol = {c: self.getValue(c) for c in self.getCompartmentIds()}
        V = np.array([vol.get(info['comp_of'].get(s), 1.0) for s in fids])
        rows = []
        if dep:  # moiety total T moves the dependent species: d conc_dep / d T = inv(Gamma[:, dep]) / V_dep
            G = self.getConservationMatrix()
            gcols = list(G.colnames)
            C = np.linalg.inv(np.asarray(G)[:, [gcols.index(d) for d in dep]]) / V[[fids.index(d) for d in dep]][:, None]
        for p in params:
            if p in zero:
                rows.append(np.zeros((len(times), S.shape[2])))
            elif p in cs:
                rows.append(np.einsum('tdx,d->tx', S[:, len(rp):], C[:, cs.index(p)]))
            else:
                rows.append(S[:, rp.index(p)])
        S = np.stack(rows, axis=1) if rows else S[:, :0]
        if k:
            S[0] = 0.0  # as libroadrunner: a zero first row for k = 1 too
        for s, x in zip(states, X[-1]):  # the model is left at the end
            if s in fids or s in rr_ids:
                self.setValue('[%s]' % s if s in fids else s, x)
        self.setValue('time', times[-1])
        idx = [states.index(s) for s in species]
        Vs = np.array([V[fids.index(s)] if s in fids else 1.0 for s in species])
        self._ran('timeSeriesSensitivities')
        return times, S[:, :, idx] * Vs, params, species  # amount sensitivities, as libroadrunner

    # ------------------------------------------------------------ reset, state
    def reset(self, options=None):
        """reset(): floating species, rate-rule values and time to their initial values;
        reset(SelectionRecord flags): also parameters (_GLOBAL_PARAMETER) or everything (ALL)."""
        if options is None:
            return self._rr.reset()
        if options == SelectionRecord.ALL:
            return self._rr.resetAll()
        if options & SelectionRecord._GLOBAL_PARAMETER:
            self._rr.resetParameter()
        if options & (SelectionRecord.TIME | SelectionRecord.RATE | SelectionRecord.FLOATING):
            self._rr.reset()

    def resetAll(self):
        self._rr.resetAll()

    def resetToOrigin(self):
        self._rr.resetToOrigin()
        self._changed()

    def resetParameter(self):
        self._rr.resetParameter()

    def saveState(self, filename, opt='b'):
        """Save the model, its values, time and event state to a file (``opt`` is ignored: one format)."""
        self._rr.saveState(filename)

    def loadState(self, filename):
        self._rr.loadState(filename)
        self._changed()

    def saveStateS(self, opt='b'):
        """saveState to a string."""
        fd, path = tempfile.mkstemp()
        os.close(fd)
        try:
            self._rr.saveState(path)
            with open(path, encoding='utf-8') as f:
                return f.read()
        finally:
            os.remove(path)

    def loadStateS(self, state):
        fd, path = tempfile.mkstemp()
        with os.fdopen(fd, 'w', encoding='utf-8') as f:
            f.write(state)
        try:
            self.loadState(path)
        finally:
            os.remove(path)

    def __getstate__(self):  # pickle / copy.deepcopy: model state and settings
        ig, ss = self.integrator, self.steadyStateSolver
        return dict(state=self.saveStateS() if self.isModelLoaded() else None,
                    integrators={n: Solver(self._rr, False, n).getSettingsMap() for n in INTEGRATORS},
                    integrator=ig.getName(), ss=ss.getName(), ss_settings=ss.getSettingsMap(),
                    sel=self.timeCourseSelections if self.isModelLoaded() else None,
                    ssel=self.steadyStateSelections if self.isModelLoaded() else None,
                    diff=self.getDiffStepSize(), flux=self.getFluxThreshold(), sens=dict(self._sens._s))

    def __setstate__(self, d):
        self.__init__()
        if d['state'] is not None:
            self.loadStateS(d['state'])
            self.timeCourseSelections, self.steadyStateSelections = d['sel'], d['ssel']
        for n, s in d['integrators'].items():
            for k, v in s.items():
                Solver(self._rr, False, n).setValue(k, v)
        self.setIntegrator(d['integrator'])
        self.setSteadyStateSolver(d['ss'])
        for k, v in d['ss_settings'].items():
            self.steadyStateSolver.setValue(k, v)
        self.setDiffStepSize(d['diff'])
        self.setFluxThreshold(d['flux'])
        self._sens._s.update(d['sens'])

    # ------------------------------------------------------------ steady state, MCA
    def steadyState(self, dict=None):
        """Solve for the steady state (``dict``: steady-state solver settings for this call);
        returns the norm of the rates of change there."""
        ss = self.steadyStateSolver
        old = {k: ss.getValue(k) for k in (dict or {})}
        for k, v in (dict or {}).items():
            ss.setValue(k, v)
        try:
            v = self._rr.steadyState()
        finally:
            for k, x in old.items():
                ss.setValue(k, x)
        self._ran('steadyState')
        return v

    def mcaSteadyState(self, *args):
        return self.steadyState()

    def getSteadyStateValues(self):
        v = np.array(self._rr.getSteadyStateValues())
        self._ran('steadyState')
        return v

    def getSteadyStateValuesNamedArray(self):
        return NamedArray(self.getSteadyStateValues()[None], colnames=self.steadyStateSelections)

    def steadyStateNamedArray(self, dict=None):
        self.steadyState(dict)
        return NamedArray(np.array(self._rr.getValues(self.steadyStateSelections))[None], colnames=self.steadyStateSelections)

    def _mat(self, method, *args):
        out = _named(getattr(self._rr, method)(*args))
        self._ran(method)
        return out

    def getFullStoichiometryMatrix(self):
        return self._mat('getFullStoichiometryMatrix')

    def getReducedStoichiometryMatrix(self):
        return self._mat('getNrMatrix')

    getNrMatrix = getReducedStoichiometryMatrix

    def getLinkMatrix(self):
        return self._mat('getLinkMatrix')

    def getL0Matrix(self):
        return self._mat('getL0Matrix')

    def getConservationMatrix(self):
        return self._mat('getConservationMatrix')

    def getFullJacobian(self):
        return self._mat('getFullJacobian')

    def getReducedJacobian(self, h=None):
        return self._mat('getReducedJacobian')

    def getFullEigenValues(self):
        self._ran('getFullEigenValues')
        return np.asarray(self._rr.getFullEigenValues())

    def getReducedEigenValues(self):
        self._ran('getReducedEigenValues')
        return np.asarray(self._rr.getReducedEigenValues())

    def getFullEigenValuesNamedArray(self):
        e = self.getFullEigenValues()
        return NamedArray(np.column_stack([e.real, e.imag]), self._rr.getFullJacobian()[1], ['real', 'imaginary'])

    def getReducedEigenValuesNamedArray(self):
        e = self.getReducedEigenValues()
        return NamedArray(np.column_stack([e.real, e.imag]), self._rr.getReducedJacobian()[1], ['real', 'imaginary'])

    def getUnscaledElasticityMatrix(self):
        return self._mat('getElasticityMatrix', False)

    def getScaledElasticityMatrix(self):
        return self._mat('getElasticityMatrix', True)

    def getUnscaledConcentrationControlCoefficientMatrix(self):
        return self._mat('getConcentrationControlCoefficientMatrix', False)

    def getScaledConcentrationControlCoefficientMatrix(self):
        return self._mat('getConcentrationControlCoefficientMatrix', True)

    def getUnscaledFluxControlCoefficientMatrix(self):
        return self._mat('getFluxControlCoefficientMatrix', False)

    def getScaledFluxControlCoefficientMatrix(self):
        return self._mat('getFluxControlCoefficientMatrix', True)

    def getCC(self, variable, parameter):
        self._ran('getCC')
        return self._rr.getCC(variable, parameter)

    def getuCC(self, variable, parameter):
        self._ran('getuCC')
        return self._rr.getuCC(variable, parameter)

    def getEE(self, reactionName, parameterName, computeSteadyState=True):
        return self._rr.getEE(reactionName, parameterName, bool(computeSteadyState))

    def getuEE(self, reactionName, parameterName, computeSteadyState=True):
        return self._rr.getuEE(reactionName, parameterName, bool(computeSteadyState))

    def getScaledFloatingSpeciesElasticity(self, reactionId, speciesId):
        return self._rr.getEE(reactionId, speciesId, False)

    def getUnscaledSpeciesElasticity(self, reactionIndex, speciesIndex):
        return self._rr.getuEE(self._rr.getReactionIds()[reactionIndex], self._rr.getFloatingSpeciesIds()[speciesIndex], False)

    def getUnscaledParameterElasticity(self, reactionName, parameterName):
        return self._rr.getuEE(reactionName, parameterName, False)

    # ------------------------------------------------------------ ids and value vectors
    def getIndependentFloatingSpeciesIds(self):
        return list(self._rr.getIndependentFloatingSpeciesIds())

    def getDependentFloatingSpeciesIds(self):
        return list(self._rr.getDependentFloatingSpeciesIds())

    def getNumIndFloatingSpecies(self):
        return len(self._rr.getIndependentFloatingSpeciesIds())

    def getNumDepFloatingSpecies(self):
        return len(self._rr.getDependentFloatingSpeciesIds())

    getNumberOfIndependentSpecies, getNumberOfDependentSpecies = getNumIndFloatingSpecies, getNumDepFloatingSpecies

    def getConservedMoietyIds(self):
        return list(self._rr.getConservedMoietyIds())

    def getConservedMoietyValues(self):
        return self.model.getConservedMoietyValues()

    def getEigenValueIds(self):
        return list(self._rr.getEigenValueIds())

    def getRatesOfChangeIds(self):
        return list(self._rr.getRatesOfChangeIds())

    def getRateRuleIds(self):
        return list(self._rr.getRateRuleIds())

    def getAssignmentRuleIds(self):
        return self.model.getAssignmentRuleIds()

    def getInitialAssignmentIds(self):
        return self.model.getInitialAssignmentIds()

    def getFloatingSpeciesInitialConcentrationIds(self):
        return self.model.getFloatingSpeciesInitConcentrationIds()

    def getAvailableTimeCourseSymbols(self):
        return self._rr.getAvailableTimeCourseSymbols()

    def getRatesOfChange(self):
        return np.array(self._rr.getRatesOfChange())

    def getRatesOfChangeNamedArray(self):
        return NamedArray(self.getRatesOfChange()[None], colnames=self.getRatesOfChangeIds())

    def getKineticLaw(self, rid):
        return roadrunner_nitro.sbml_edit('get_kinetic_law', self.getSBML(), rid)

    # ------------------------------------------------------------ model editing
    def _edit(self, op, args, regen=True):
        self._set_sbml(roadrunner_nitro.sbml_edit(op, self._rr.getSBML(), *args), regen)

    def _set_sbml(self, sbml, regen=True):
        self._rr.setSBML(sbml, bool(regen))
        self._changed()

    def regenerateModel(self, forceRegenerate=True, reset=False):
        """Recompile the (edited) model, keeping the current values (``reset``: then resetAll)."""
        self._rr.regenerate()
        self._changed()
        if reset:
            self.resetAll()

    def addSpecies(self, sid, compartment, initAmount=0.0, initConcentration=0.0, hasOnlySubstanceUnits=False,
                   boundaryCondition=False, substanceUnits='', forceRegenerate=True):
        if initConcentration and not initAmount:
            self.addSpeciesConcentration(sid, compartment, initConcentration, hasOnlySubstanceUnits, boundaryCondition,
                                         substanceUnits, forceRegenerate)
        else:
            self.addSpeciesAmount(sid, compartment, initAmount, hasOnlySubstanceUnits, boundaryCondition,
                                  substanceUnits, forceRegenerate)

    def addReaction(self, *args, **kwargs):
        """addReaction(sbmlRep[, forceRegenerate]) or addReaction(rid, reactants, products, kineticLaw[, forceRegenerate])."""
        regen = kwargs.pop('forceRegenerate', True)
        if len(args) in (1, 2) and not kwargs:
            if len(args) == 2:
                regen = args[1]
            return self._edit('add_reaction_from_sbml', [args[0]], regen)
        a = dict(zip(['rid', 'reactants', 'products', 'kineticLaw', 'forceRegenerate'], args))
        a.update(kwargs)
        self._edit('add_reaction', [a['rid'], list(a['reactants']), list(a['products']), a['kineticLaw']],
                   a.get('forceRegenerate', regen))

    def removeRules(self, vid, useInitialValue=False, forceRegenerate=True):
        """Remove the rules of ``vid`` (``useInitialValue``: its value becomes the rule's value at time 0)."""
        v = self.getValue('init(%s)' % vid) if useInitialValue else None
        self._edit('remove_rules', [vid], forceRegenerate and v is None)
        if v is not None:
            self.setInitValue(vid, v, forceRegenerate)

    def setInitAmount(self, sid, initAmount, forceRegenerate=True):
        self._edit('set_init_amount', [sid, float(initAmount)], forceRegenerate)

    def setInitConcentration(self, sid, initConcentration, forceRegenerate=True):
        self._edit('set_init_concentration', [sid, float(initConcentration)], forceRegenerate)

    def setInitValue(self, sid, initValue, forceRegenerate=True):
        """Initial value of any species (amount), compartment or parameter."""
        sbml = self._rr.getSBML()
        if sid in self.getFloatingSpeciesIds() + self.getBoundarySpeciesIds():
            return self.setInitAmount(sid, initValue, forceRegenerate)
        tag, attr = ('compartment', 'size') if sid in self.getCompartmentIds() else ('parameter', 'value')
        self._set_sbml(_set_attr(sbml, tag, sid, attr, repr(float(initValue))), forceRegenerate)

    def removeInitialAssignment(self, vid, forceRegenerate=True, errIfNotExist=True, replaceInitVal=True):
        """Remove the initial assignment of ``vid`` (``replaceInitVal``: its value becomes the initial value).
        Returns whether the model was regenerated."""
        sbml = self._rr.getSBML()
        pat = re.compile(r'<(?:\w+:)?initialAssignment\b[^>]*?\bsymbol\s*=\s*"%s".*?</(?:\w+:)?initialAssignment>\s*'
                         % re.escape(vid), re.S)
        if not pat.search(sbml):
            if errIfNotExist:
                raise ValueError("no initial assignment for '%s'" % vid)
            return False
        v = self.getValue('init(%s)' % vid) if replaceInitVal else None
        self._set_sbml(pat.sub('', sbml, count=1), forceRegenerate and v is None)
        if v is not None:
            self.setInitValue(vid, v, forceRegenerate)
        return bool(forceRegenerate)

    def setInitStoichiometryValue(self, speciesId, reactionId, initValue):
        """Stoichiometry of the (unnamed) species reference of ``speciesId`` in reaction ``reactionId``."""
        sbml = self._rr.getSBML()
        m = re.search(r'<(?:\w+:)?reaction\b[^>]*?\bid\s*=\s*"%s".*?</(?:\w+:)?reaction>' % re.escape(reactionId), sbml, re.S)
        if m is None:
            raise ValueError("no reaction '%s'" % reactionId)
        block = _set_attr(m.group(0), 'speciesReference', speciesId, 'stoichiometry', repr(float(initValue)), 'species')
        self._set_sbml(sbml[:m.start()] + block + sbml[m.end():])

    def setModelId(self, id):
        self._set_sbml(_set_attr(self._rr.getSBML(), 'model', None, 'id', id), False)

    def setModelName(self, name):
        self._set_sbml(_set_attr(self._rr.getSBML(), 'model', None, 'name', name), False)

    def setGlobalParameterByName(self, param, value):
        self.setValue(param, value)

    def getGlobalParameterByName(self, param):
        return self.getValue(param)

    # The model getters libroadrunner's tellurium attaches to RoadRunner (r.getFloatingSpeciesIds() ...)
    _model_functions = [
        'getBoundarySpeciesIds', 'getBoundarySpeciesConcentrationIds', 'getBoundarySpeciesConcentrations',
        'getBoundarySpeciesAmounts', 'getNumBoundarySpecies',
        'getFloatingSpeciesIds', 'getFloatingSpeciesConcentrationIds', 'getFloatingSpeciesConcentrations',
        'getFloatingSpeciesAmounts', 'getNumFloatingSpecies',
        'getGlobalParameterIds', 'getGlobalParameterValues', 'getNumGlobalParameters',
        'getCompartmentIds', 'getCompartmentVolumes', 'getNumCompartments',
        'getNumConservedMoieties', 'getReactionIds', 'getReactionRates', 'getNumReactions', 'getNumEvents',
        'getNumRateRules']

    # ---------------------------------------------------------------------
    # Jarnac compatibility layer
    # ---------------------------------------------------------------------
    def fjac(self):
        return self.getFullJacobian()

    def sm(self):
        return self.getFullStoichiometryMatrix()

    def rs(self):
        return self.model.getReactionIds()

    def fs(self):
        return self.model.getFloatingSpeciesIds()

    def bs(self):
        return self.model.getBoundarySpeciesIds()

    def ps(self):
        return self.model.getGlobalParameterIds()

    def vs(self):
        return self.model.getCompartmentIds()

    def dv(self):
        return self.getRatesOfChange()

    def rv(self):
        return self.model.getReactionRates()

    def sv(self):
        return self.model.getFloatingSpeciesConcentrations()

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

        from tellurium_nitro import SBMLDiagram
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


for _n, (_op, _k, _d) in _EDITS.items():
    setattr(ExtendedRoadRunner, _n, _edit_method(_n, _op, _k, _d))


def _model_function(key):
    def f(self):
        return getattr(self.model, key)()
    f.__name__ = key
    return f


for _n in ExtendedRoadRunner._model_functions:
    setattr(ExtendedRoadRunner, _n, _model_function(_n))


def _vector(ids, fmt, named):  # r.get<Independent|Dependent|>FloatingSpecies<Amounts|Concentrations><V|NamedArray>()
    def f(self):
        names = getattr(self._rr, ids)()
        v = np.array(self._rr.getValues([fmt % s for s in names]))
        return NamedArray(v[None], colnames=names) if named else v
    return f


for _g, _ids in (('', 'getFloatingSpeciesIds'), ('Independent', 'getIndependentFloatingSpeciesIds'),
                 ('Dependent', 'getDependentFloatingSpeciesIds'), ('Boundary', 'getBoundarySpeciesIds')):
    _sp = 'BoundarySpecies' if _g == 'Boundary' else _g + 'FloatingSpecies'
    for _q, _fmt in (('Amounts', '%s'), ('Concentrations', '[%s]')):
        for _sfx, _named_ in (('V', False), ('NamedArray', True)):
            if _g in ('', 'Boundary') and _sfx == 'V':
                continue
            _f = _vector(_ids, _fmt, _named_)
            _f.__name__ = 'get%s%s%s' % (_sp, _q, _sfx)
            setattr(ExtendedRoadRunner, _f.__name__, _f)
    if _g in ('Independent', 'Dependent'):
        _f = _vector(_ids, "%s'", False)
        _f.__name__ = 'get%sRatesOfChange' % _g
        setattr(ExtendedRoadRunner, _f.__name__, _f)
        _f = _vector(_ids, "%s'", True)
        _f.__name__ = 'get%sRatesOfChangeNamedArray' % _g
        setattr(ExtendedRoadRunner, _f.__name__, _f)

RoadRunner = ExtendedRoadRunner
