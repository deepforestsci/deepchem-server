"""Text to AST.

Lark owns tokenisation and indentation; this module owns the translation into
our own AST so no Lark object reaches the validator or interpreter.
"""

import functools
import json
import os
from typing import Any, List

from lark import Lark, Token, Tree
from lark.exceptions import LarkError
from lark.indenter import Indenter

from deepchem_server.core.script import ast_nodes as ast
from deepchem_server.core.script.errors import ScriptParseError


GRAMMAR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "grammar.lark")


class ScriptIndenter(Indenter):
    NL_type = "_NL"
    OPEN_PAREN_types = ["LPAR", "LSQB", "LBRACE"]
    CLOSE_PAREN_types = ["RPAR", "RSQB", "RBRACE"]
    INDENT_type = "_INDENT"
    DEDENT_type = "_DEDENT"
    tab_len = 8


@functools.lru_cache(maxsize=None)
def _parser() -> Lark:
    with open(GRAMMAR_PATH) as handle:
        grammar = handle.read()
    # maybe_placeholders defaults to True in Lark >=1.0, which makes an empty
    # "[...]" optional group (used by call/list_literal/dict_literal for their
    # empty-body case) yield a single ``None`` child instead of no children at
    # all. Disable it so an empty group contributes zero children, matching
    # what every builder below assumes.
    return Lark(
        grammar,
        parser="lalr",
        postlex=ScriptIndenter(),
        propagate_positions=True,
        maybe_placeholders=False,
    )


_MATCHING_OPENER = {")": "(", "]": "[", "}": "{"}
_OPENERS = "".join(_MATCHING_OPENER.values())
_CLOSERS = "".join(_MATCHING_OPENER)


def _find_unclosed_opener(text: str) -> Any:
    """Return ``(char, line, column)`` of the earliest bracket never closed anywhere in ``text``, else None.

    LALR reports an unclosed bracket where parsing diverges, often a line
    later. Scanning the whole script (skipping string literals) only blames
    an opener that genuinely never closes; a bracket-kind mismatch defers to
    Lark's own message.
    """
    stack: List[Any] = []
    in_string = False
    escape = False
    line, column = 1, 0
    for ch in text:
        column += 1
        if in_string:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_string = False
        elif ch == '"':
            in_string = True
        elif ch in _OPENERS:
            stack.append((ch, line, column))
        elif ch in _CLOSERS:
            if stack:
                if stack[-1][0] == _MATCHING_OPENER[ch]:
                    stack.pop()
                else:
                    return None
        if ch == "\n":
            line += 1
            column = 0
            in_string = False
            escape = False
    return stack[0] if stack else None


def parse_script(text: str) -> ast.Script:
    """Parse script text into a Script, raising ScriptParseError with a position."""
    if not text.endswith("\n"):
        text = text + "\n"
    try:
        tree = _parser().parse(text)
    except LarkError as exc:
        opener = _find_unclosed_opener(text)
        if opener is not None:
            char, line, column = opener
            raise ScriptParseError(f"{char!r} was never closed", line, column) from exc
        line = getattr(exc, "line", 0) or 0
        column = getattr(exc, "column", 0) or 0
        raise ScriptParseError(str(exc).strip().split("\n")[0], line, column) from exc
    return ast.Script(statements=[_statement(child) for child in tree.children])


def _pos(node: Any) -> dict:
    """Return the {line, column} a human would point at for ``node``.

    ``node`` is either a Lark ``Token`` (which carries ``.line``/``.column``
    natively) or a Lark ``Tree`` (whose position lives on ``.meta``, and only
    when Lark could determine one -- ``meta.empty`` is True otherwise). For an
    unpositioned Tree, fall back to the first descendant Token's position.
    """
    if isinstance(node, Tree):
        meta = node.meta
        if not meta.empty:
            return {"line": meta.line or 0, "column": meta.column or 0}
        token = _first_token(node)
        if token is not None:
            return {"line": token.line or 0, "column": token.column or 0}
        return {"line": 0, "column": 0}
    return {"line": getattr(node, "line", 0) or 0, "column": getattr(node, "column", 0) or 0}


def _first_token(node: Tree) -> Any:
    """Depth-first search for the first Token descendant of a Tree, or None."""
    for child in node.children:
        if isinstance(child, Token):
            return child
        if isinstance(child, Tree):
            token = _first_token(child)
            if token is not None:
                return token
    return None


def _statement(node: Tree) -> Any:
    handler = _STATEMENTS.get(node.data)
    if handler is None:
        raise ScriptParseError(f"unsupported statement {node.data!r}", **_pos(node))
    return handler(node)


def _assign(node: Tree) -> ast.Assign:
    name, value = node.children
    return ast.Assign(name=str(name), value=_expr(value), **_pos(name))


def _call_stmt(node: Tree) -> ast.CallStmt:
    call = _expr(node.children[0])
    return ast.CallStmt(call=call, line=call.line, column=call.column)


def _rank_on(node: Tree) -> ast.RankOn:
    value = _expr(node.children[0])
    return ast.RankOn(value=value, line=value.line, column=value.column)


def _body(block: Tree) -> List[Any]:
    return [_statement(child) for child in block.children]


def _if_stmt(node: Tree) -> ast.If:
    branches = [ast.Branch(test=_expr(node.children[0]), body=_body(node.children[1]))]
    orelse: List[Any] = []
    for child in node.children[2:]:
        if child.data == "elif_clause":
            branches.append(ast.Branch(test=_expr(child.children[0]), body=_body(child.children[1])))
        else:
            orelse = _body(child.children[0])
    return ast.If(branches=branches, orelse=orelse, **_pos(node.children[0]))


def _for_stmt(node: Tree) -> ast.For:
    target, iterable, block = node.children
    return ast.For(target=str(target), iterable=_expr(iterable), body=_body(block), **_pos(target))


def _expr(node: Any) -> Any:
    if isinstance(node, Token):
        return _token(node)
    builder = _EXPRESSIONS.get(node.data)
    if builder is None:
        raise ScriptParseError(f"unsupported expression {node.data!r}", **_pos(node))
    return builder(node)


def _decode_string(token: Any, what: str) -> str:
    """Decode a STRING token (quotes included); an escape JSON rejects is a positioned parse error."""
    try:
        return json.loads(str(token))
    except json.JSONDecodeError as exc:
        raise ScriptParseError(f"invalid escape in {what}: {exc.msg}", **_pos(token)) from exc


_KEYWORDS = {"TRUE": True, "FALSE": False, "NONE": None}


def _token(token: Token) -> Any:
    position = _pos(token)
    if token.type == "STRING":
        return ast.Literal(value=_decode_string(token, "string literal"), **position)
    if token.type == "SIGNED_NUMBER":
        text = str(token)
        return ast.Literal(value=float(text) if "." in text else int(text), **position)
    if token.type in _KEYWORDS:
        return ast.Literal(value=_KEYWORDS[token.type], **position)
    raise ScriptParseError(f"unsupported token {token.type}", **position)


def _dict(node: Tree) -> ast.DictExpr:
    items = {}
    for pair in node.children:
        key, value = pair.children
        items[_decode_string(key, "dict key")] = _expr(value)
    return ast.DictExpr(items=items, **_pos(node))


def _call(node: Tree) -> ast.Call:
    name = node.children[0]
    kwargs = {}
    for pair in node.children[1:]:
        kwargs[str(pair.children[0])] = _expr(pair.children[1])
    return ast.Call(name=str(name), kwargs=kwargs, **_pos(name))


def _atom_expr(node: Tree) -> Any:
    value = _expr(node.children[0])
    for trailer in node.children[1:]:
        if trailer.data == "attr_trailer":
            value = ast.Attribute(value=value, attr=str(trailer.children[0]), line=value.line, column=value.column)
        else:
            value = ast.Index(value=value, index=_expr(trailer.children[0]), line=value.line, column=value.column)
    return value


def _binary(node: Tree) -> Any:
    children = node.children
    builder = ast.Compare if node.data == "comparison" else ast.BinOp
    value = _expr(children[0])
    for op, right in zip(children[1::2], children[2::2]):
        value = builder(op=str(op), left=value, right=_expr(right), line=value.line, column=value.column)
    return value


_STATEMENTS = {
    "assign": _assign,
    "call_stmt": _call_stmt,
    "if_stmt": _if_stmt,
    "for_stmt": _for_stmt,
    "rank_on_stmt": _rank_on,
}

_EXPRESSIONS = {
    "call": _call,
    "name_ref": lambda n: ast.Name(id=str(n.children[0]), **_pos(n.children[0])),
    "list_literal": lambda n: ast.ListExpr(items=[_expr(c) for c in n.children], **_pos(n)),
    "dict_literal": _dict,
    "range_form": lambda n: ast.RangeForm(count=int(n.children[0]), **_pos(n)),
    "atom_expr": _atom_expr,
    "not_op": lambda n: ast.UnaryOp(op="not", operand=_expr(n.children[0]), **_pos(n)),
    "neg": lambda n: ast.UnaryOp(op="-", operand=_expr(n.children[0]), **_pos(n)),
    "comparison": _binary,
    "arith": _binary,
    "term": _binary,
    "or_test": lambda n: ast.BoolOp(op="or", values=[_expr(c) for c in n.children], **_pos(n)),
    "and_test": lambda n: ast.BoolOp(op="and", values=[_expr(c) for c in n.children], **_pos(n)),
}
