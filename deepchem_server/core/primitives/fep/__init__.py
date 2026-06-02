# flake8: noqa
import warnings

from deepchem_server.core.primitives.fep.rbfe.collate_rbfe_results import collate_rbfe_results
from deepchem_server.core.primitives.fep.rbfe.run_rbfe import run_rbfe

# Suppress noisy third-party deprecation warnings from openff/gufe dependencies
# without replacing warnings.warn globally (which would break catch_warnings).
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"gufe\..*")
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"openff\..*")
