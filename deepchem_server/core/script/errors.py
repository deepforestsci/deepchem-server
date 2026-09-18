"""Errors raised while parsing, validating or executing a workflow script.

Every error carries a source position because the Executor agent repairs
scripts from these messages; an error without a line number is not
actionable.
"""


class ScriptError(Exception):
    """Base class for workflow script failures."""

    kind = "script"

    def __init__(self, message: str, line: int = 0, column: int = 0):
        super().__init__(message)
        self.message = message
        self.line = line
        self.column = column

    def to_dict(self) -> dict:
        return {"kind": self.kind, "line": self.line, "column": self.column, "message": self.message}


class ScriptParseError(ScriptError):
    """The text is not a valid script."""

    kind = "parse"


class ScriptValidationError(ScriptError):
    """The script parses but cannot be executed as written."""

    kind = "validation"


class ScriptRuntimeError(ScriptError):
    """Execution began and could not complete."""

    kind = "runtime"
