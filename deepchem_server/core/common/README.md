# DeepChem Core Common

Common utilities for DeepChem Server including:

- **Address**: DeepChem address parsing and manipulation (`deepchem://profile/project/key`)
- **Cards**: DataCard and ModelCard for metadata management
- **Config**: Global configuration management (datastore settings)
- **DataStore**: Abstract and disk-based data storage implementations
- **Progress Logger**: Logging utilities for job progress tracking
- **Utils**: Helper functions for job execution and data parsing

## Installation

### Development (UV Workspace)

When working within the monorepo, use UV workspace:

```bash
# From repository root
uv sync
```

### Standalone Installation

```bash
pip install deepchem-core-common
```

### With dev dependencies

```bash
pip install deepchem-core-common[dev]
```

## Usage

```python
from deepchem_server.core.common import (
    DeepchemAddress,
    DataCard,
    ModelCard,
    DiskDataStore,
    get_datastore,
    set_datastore,
)

# Parse a DeepChem address
address = DeepchemAddress("deepchem://profile/project/data.csv")
print(address.profile, address.project, address.key)

# Create a DataCard
card = DataCard(
    address="",
    file_type="csv",
    data_type="pandas.DataFrame",
    description="My dataset"
)
```

## Development

```bash
# Install with dev dependencies
uv pip install -e ".[dev]"

# Run tests
pytest
```

