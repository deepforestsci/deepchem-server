"""AST for the workflow scripting language.

Nodes are plain dataclasses rather than a class hierarchy with behaviour so
the validator and interpreter stay separable: neither can smuggle execution
into the tree itself.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union


@dataclass
class Node:
    line: int = 0
    column: int = 0


@dataclass
class Literal(Node):
    value: Any = None


@dataclass
class Name(Node):
    id: str = ""


@dataclass
class ListExpr(Node):
    items: List[Any] = field(default_factory=list)


@dataclass
class DictExpr(Node):
    items: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Attribute(Node):
    value: Any = None
    attr: str = ""


@dataclass
class Index(Node):
    value: Any = None
    index: Any = None


@dataclass
class Call(Node):
    name: str = ""
    kwargs: Dict[str, Any] = field(default_factory=dict)


@dataclass
class UnaryOp(Node):
    op: str = ""
    operand: Any = None


@dataclass
class BinOp(Node):
    op: str = ""
    left: Any = None
    right: Any = None


@dataclass
class Compare(Node):
    op: str = ""
    left: Any = None
    right: Any = None


@dataclass
class BoolOp(Node):
    op: str = ""
    values: List[Any] = field(default_factory=list)


@dataclass
class RangeForm(Node):
    count: int = 0


@dataclass
class Assign(Node):
    name: str = ""
    value: Any = None


@dataclass
class CallStmt(Node):
    call: Optional[Call] = None


@dataclass
class Branch:
    test: Any
    body: List[Any]


@dataclass
class If(Node):
    branches: List[Branch] = field(default_factory=list)
    orelse: List[Any] = field(default_factory=list)


@dataclass
class For(Node):
    target: str = ""
    iterable: Any = None
    body: List[Any] = field(default_factory=list)


@dataclass
class RankOn(Node):
    value: Any = None


Statement = Union[Assign, CallStmt, If, For, RankOn]


@dataclass
class Script:
    statements: List[Statement] = field(default_factory=list)
