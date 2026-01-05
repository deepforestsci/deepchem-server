"""
DEPRECATED: This module is deprecated. Import from specific packages instead.

- API utilities: from deepchem_server.api.utils import ...
- Job utilities: from deepchem_server.services.jobs import run_job

This module is kept for backward compatibility and will emit deprecation warnings.
"""
import warnings


warnings.warn(
    "deepchem_server.utils is deprecated. "
    "Use deepchem_server.api.utils or deepchem_server.services.jobs instead.",
    DeprecationWarning,
    stacklevel=2,
)

# Re-exports for backward compatibility
from deepchem_server.api.utils import (  # noqa: E402
    parse_boolean_none_values_from_kwargs, parse_dict_with_datatypes,
)
from deepchem_server.api.utils import (  # noqa: E402
    upload_data as _upload_data,)
from deepchem_server.services.jobs import run_job  # noqa: E402
from deepchem_server.services.jobs.utils import (  # noqa: E402
    init_datastore as _init_datastore,)
from deepchem_server.services.jobs.utils import run_job_sync  # noqa: E402


__all__ = [
    "_upload_data",
    "_init_datastore",
    "run_job",
    "run_job_sync",
    "parse_boolean_none_values_from_kwargs",
    "parse_dict_with_datatypes",
]
