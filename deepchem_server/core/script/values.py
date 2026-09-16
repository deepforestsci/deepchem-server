"""Runtime values a script can hold.

A primitive result is deliberately not a dict: exposing a fixed set of named
fields means a script has no route to interpreter internals, and the failure
modes of `.address` are explicit rather than an IndexError.
"""

from dataclasses import dataclass, field
from typing import Dict, List

from deepchem_server.core.script.errors import ScriptRuntimeError


@dataclass
class ResultValue:
    """The outcome of one primitive call, as seen from inside a script."""

    program_name: str
    addresses: List[str] = field(default_factory=list)
    metrics: Dict[str, float] = field(default_factory=dict)
    kind: str = "address"

    @property
    def address(self) -> str:
        """The single produced address.

        Raises when the call produced none or several, because silently
        picking one would make a script's meaning depend on a primitive's
        arity rather than on what the script says.
        """
        if not self.addresses:
            raise ScriptRuntimeError(f"{self.program_name} produced no address")
        if len(self.addresses) > 1:
            raise ScriptRuntimeError(f"{self.program_name} produced {len(self.addresses)} addresses; "
                                     "index one explicitly with .addresses[i]")
        return self.addresses[0]


@dataclass(frozen=True)
class Limits:
    """Ceilings the server enforces. A script cannot read or change these."""

    max_primitive_calls: int = 64
    max_total_iterations: int = 256
    max_loop_iterations: int = 64
    max_nesting_depth: int = 4
    max_statements: int = 500
    max_wall_time_s: float = 3600.0
    max_script_bytes: int = 65536


DEFAULT_LIMITS = Limits()
