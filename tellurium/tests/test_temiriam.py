"""
Testing temiriam module
"""
from __future__ import absolute_import, print_function
import pytest
import requests
import roadrunner
from tellurium import temiriam


def _fetch(urn):
    try:
        return temiriam.getSBMLFromBiomodelsURN(urn)
    except requests.exceptions.RequestException as e:
        pytest.skip('BioModels not reachable (network test): {}'.format(e))


def test_getSBMLFromBiomodelsURN1():
    """ Check that string is returned.

    :return:
    """
    urn = 'urn:miriam:biomodels.db:BIOMD0000000139'
    sbml = _fetch(urn)
    assert sbml is not None
    # check that string
    assert isinstance(sbml, str)


def test_getSBMLFromBiomodelsURN2():
    """ Check that model can be loaded in roadrunner.

    :return:
    """
    urn = 'urn:miriam:biomodels.db:BIOMD0000000139'
    sbml = _fetch(urn)

    print("*" * 80)
    print(type(sbml))
    print("*" * 80)
    print(sbml)
    print("*" * 80)

    r = roadrunner.RoadRunner(sbml)
    assert r is not None
