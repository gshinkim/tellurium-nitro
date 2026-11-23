"""
Helpers to work with SBO terms in antimony.
"""

from __future__ import print_function, division, absolute_import

import re
import xml.etree.ElementTree as ET

from .antimony_regex import getModelStartRegex, getModelEndRegex, getFunctionStartRegex, getSBORegex, getFunctionSBORegex

def filterIfEmpty(l):
    """ If l is one line (comment), filter."""
    if len(l) == 1:
        return []
    else:
        return l

class SBOError(RuntimeError):
    pass

class antimonySBOConverter(object):
    def __init__(self, doc):
        """ Create an SBO converter.

        :param doc: the SBML document's root element (xml.etree.ElementTree)
        """
        self.doc = doc
        self.model = _child(doc, 'model')
        if self.model is None:
            raise RuntimeError('No SBML model')

    def convert(self, antimony_str):
        """ Add SBO terms to the Antimony string corresponding to this SBML document.

        :param antimony_str: Antimony string. Should represent the same SBML this object was initialized with.
        :return: The Antimony string with the SBO terms added.
        """

        model_start = re.compile(getModelStartRegex())
        n_model_starts = 0

        function_start = re.compile(getFunctionStartRegex())

        model_end = re.compile(getModelEndRegex())
        n_model_ends = 0
        model_end_index = None

        scopes = []

        n_leading_spaces = 0

        lines = antimony_str.splitlines()


        for n,line in enumerate(lines):
            if model_start.match(line) != None:
                n_model_starts += 1
                scopes.append('model')
            if function_start.match(line) != None:
                scopes.append('function')
            elif model_end.match(line) != None:
                if n_model_ends > 0:
                    raise SBOError('Multiple embedded models (e.g. comp) not supported for SBO converter')
                if len(scopes) == 0:
                    raise RuntimeError('Unbalanced begin/end blocks')
                scope = scopes.pop()
                if scope == 'model':
                    n_model_ends += 1
                    model_end_index = n
                elif scope == 'function':
                    pass
                else:
                    raise RuntimeError('Unknown scope')
            else:
                # calculate leading whitespace
                ws = len(line) - len(line.lstrip(' '))
                if ws > n_leading_spaces and n_leading_spaces == 0:
                    n_leading_spaces = ws
        if n_model_starts != n_model_ends:
            raise RuntimeError('Antimony model begin/end blocks unbalanced - missing begin/end marker?')
        if n_model_starts > 1:
            raise SBOError('Multiple embedded models (e.g. comp) not supported for SBO converter')
        if n_model_ends > 1:
            raise SBOError('Multiple embedded models (e.g. comp) not supported for SBO converter')
        # add SBO terms right before the end
        lead_space = ' '*n_leading_spaces
        if model_end_index is not None:
            # get all sbo terms
            sbo_terms = self.getAllSBOTerms()
            if sbo_terms:
                lines[model_end_index:model_end_index] = ['', lead_space + '// SBO terms:'] + \
                    [lead_space + term for term in sbo_terms] + ['']
        return '\n'.join(lines)

    def getAllSBOTerms(self):
        """ Get a list of all SBO terms base on the SBML passed to the constructor. """
        comps = self.getCompartmentSBOTerms()
        species = self.getSpeciesSBOTerms()
        params = self.getParameterSBOTerms()
        rxns = self.getReactionSBOTerms()
        return comps + species + params + rxns

    def _terms(self, lst, item):
        l = _child(self.model, lst)
        return [self.createSBOTermStringForElt(e) for e in ([] if l is None else l)
                if _local(e.tag) == item and self.hasSBOTerm(e)]

    def getCompartmentSBOTerms(self):
        return filterIfEmpty(['// - Compartment SBO Terms:'] + self._terms('listOfCompartments', 'compartment'))

    def getSpeciesSBOTerms(self):
        return filterIfEmpty(['// - Species SBO Terms:'] + self._terms('listOfSpecies', 'species'))

    def getParameterSBOTerms(self):
        return filterIfEmpty(['// - Parameter SBO Terms:'] + self._terms('listOfParameters', 'parameter'))

    def getReactionSBOTerms(self):
        return filterIfEmpty(['// - Reaction SBO Terms:'] + self._terms('listOfReactions', 'reaction'))

    def hasSBOTerm(self, elt):
        return elt.get('sboTerm') is not None

    def createSBOTermStringForElt(self, elt):
        if self.hasSBOTerm(elt):
            return elt.get('id') + '.sboTerm = ' + 'SBO:{:07d};'.format(self.getSBOTermForElt(elt))
        return None

    def getSBOTermForElt(self, elt):
        if not self.hasSBOTerm(elt):
            raise RuntimeError('No SBO term set.')
        return int(elt.get('sboTerm').split(':')[-1])

    @classmethod
    def fromSBMLFile(cls, sbml_file):
        """ Construct from SBML file. """
        with open(sbml_file, encoding='utf-8') as f:
            return cls.fromSBMLString(f.read())

    @classmethod
    def fromSBMLString(cls, sbml_str):
        """ Construct from SBML string. """
        try:
            return antimonySBOConverter(ET.fromstring(sbml_str.encode('utf-8')))
        except ET.ParseError as e:
            raise RuntimeError('Errors encountered loading SBML: %s' % e)


def _local(tag):
    return tag.rsplit('}', 1)[-1]


def _child(e, name):
    return next((c for c in e if _local(c.tag) == name), None)


class antimonySBOParser(object):
    def __init__(self, antimony_str):
        """Initialize from an Antimony string."""
        self.antimony_str = antimony_str

    def elideSBOTerms(self):
        """Remove SBO terms from self.antimony_str.
        Remove SBO terms for functions. See https://github.com/sys-bio/tellurium/issues/340.

        :return: Antimony string without SBO terms."""

        self.sbo_map = {}

        model_start = re.compile(getModelStartRegex())
        n_model_starts = 0

        model_end = re.compile(getModelEndRegex())
        n_model_ends = 0
        model_end_index = None

        sbo_term = re.compile(getSBORegex())

        fct_start = re.compile(getFunctionStartRegex())
        fct_end = re.compile(getModelEndRegex())
        function_sbo = re.compile(getFunctionSBORegex())
        fct_id = ''

        n_functions = 0
        in_function = False

        n_leading_spaces = 0

        lines = self.antimony_str.splitlines()

        out_lines = []
        for n,line in enumerate(lines):
            sbo_match = sbo_term.match(line)
            fct_match = fct_start.match(line)
            fct_sbo_match = function_sbo.match(line)

            if model_start.match(line) != None:
                if n_model_starts > 0:
                    raise SBOError('Multiple embedded models (e.g. comp) not supported for SBO converter')
                n_model_starts += 1
                out_lines.append(line)
            elif model_end.match(line) != None:
                if in_function:
                    in_function = False
                    fct_id = ''
                else:
                    n_model_ends += 1
                    model_end_index = n
                out_lines.append(line)
            elif fct_match != None:
                if not in_function:
                    n_functions += 1
                    in_function = True
                    out_lines.append(line)
                    fct_id = fct_match.group(1)
                    #print('Match function id {}'.format(fct_id))
                else:
                    raise RuntimeError('Nested function: {}'.format(line))
            elif sbo_term.match(line):
                elt_id = sbo_match.group(1)
                sbo = int(sbo_match.group(3))
                self.sbo_map[elt_id] = sbo
                #print('Match SBO term {}->{}'.format(elt_id,sbo))
            elif fct_sbo_match != None:
                sbo = int(fct_sbo_match.group(2))
                #self.sbo_map[elt_id] = sbo
                #print('Match function SBO term {}->{}'.format(fct_id,sbo))
                self.sbo_map[fct_id] = sbo
            else:
                out_lines.append(line)
        if n_model_starts != n_model_ends:
            raise RuntimeError('Antimony model begin/end blocks unbalanced - missing begin/end marker?')
        if n_model_starts > 1:
            raise SBOError('Multiple embedded models (e.g. comp) not supported for SBO converter')
        if n_model_ends > 1:
            raise SBOError('Multiple embedded models (e.g. comp) not supported for SBO converter')

        return '\n'.join(out_lines)

    def addSBOsToSBML(self, sbml_str):
        """Add SBO terms to an SBML string. Must have called
        elideSBOTerms first to populate self.sbo_map."""

        out_sbml = sbml_str
        for elt_id, sbo in self.sbo_map.items():  # the sboTerm attribute of the element with that id
            pat = re.compile(r'<([\w:]+)\b([^>]*?\bid\s*=\s*"%s"[^>]*?)(/?)>' % re.escape(elt_id))
            m = pat.search(out_sbml)
            if m is not None:
                body = re.sub(r'\ssboTerm\s*=\s*"[^"]*"', '', m.group(2)) + ' sboTerm="SBO:{:07d}"'.format(sbo)
                out_sbml = out_sbml[:m.start()] + '<%s%s%s>' % (m.group(1), body, m.group(3)) + out_sbml[m.end():]
        return out_sbml
