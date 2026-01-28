from enum import Enum


EMPTY_STRING = ''

LIGAND = 'ligand'
PROTEIN = 'protein'
SOLVENT = 'solvent'

SUCCESS = 'success'
FAILURE = 'failure'

SAMPLER = 'sampler'
DEBUG = 'debug'

TOLUENE = 'toluene'
BENZENE = 'benzene'


class NetworkPlanningConstants:

    class PerturbationNetworkType(Enum):
        """Enum of Perturbation Network Types

        This is used for validating the Network Type passed to the
        get_perturbation_network() method. Add the corresponding Key-Value pair here
        when adding a new Network Type.

        """
        RADIAL = 'RADIAL'
        MINIMAL_SPANNING = 'MINIMAL_SPANNING'
        MAXIMAL = 'MAXIMAL',

    class ScorerType(Enum):
        """Enum of Scorer Types

        This is used for validating the Scorer passed to the
        get_perturbation_network() method. Add the corresponding Key-Value pair here
        when adding a new Scorer Type.
        """
        LOMAP = 'LOMAP'
