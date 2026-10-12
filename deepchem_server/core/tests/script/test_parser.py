import pytest

from deepchem_server.core.script.ast_nodes import Assign, Call, Literal, RankOn, Script
from deepchem_server.core.script.errors import ScriptParseError
from deepchem_server.core.script.parser import parse_script


def test_parses_assignment_and_rank_on():
    script = parse_script('x = "deepchem://a/b/c.csv"\nrank_on x\n')

    assert isinstance(script, Script)
    assert len(script.statements) == 2
    assign, rank = script.statements
    assert isinstance(assign, Assign)
    assert assign.name == "x"
    assert isinstance(assign.value, Literal)
    assert assign.value.value == "deepchem://a/b/c.csv"
    assert assign.line == 1
    assert isinstance(rank, RankOn)
    assert rank.line == 2


def test_parses_keyword_only_call():
    script = parse_script('f = featurize(dataset_address="a", featurizer="ecfp")\nrank_on f\n')

    call = script.statements[0].value
    assert isinstance(call, Call)
    assert call.name == "featurize"
    assert sorted(call.kwargs) == ["dataset_address", "featurizer"]
    assert call.kwargs["featurizer"].value == "ecfp"


def test_positional_argument_is_a_parse_error_with_position():
    with pytest.raises(ScriptParseError) as excinfo:
        parse_script('f = featurize("a")\nrank_on f\n')

    # The "(" closes on the same line, so this is Lark's own diagnosis, not
    # the unclosed-bracket heuristic.
    assert excinfo.value.line == 1
    assert excinfo.value.column == 15
    assert "never closed" not in excinfo.value.message
    assert "Unexpected token" in excinfo.value.message


def test_invalid_escape_in_string_literal_raises_parse_error_with_position():
    # ESCAPED_STRING tokenizes "\q" but JSON decoding rejects it.
    with pytest.raises(ScriptParseError) as excinfo:
        parse_script("x = 1\ny = 2\nz = \"\\q\"\nrank_on z\n")

    assert excinfo.value.line == 3


def test_unclosed_opener_reports_its_own_line_not_the_next():
    # Lark notices the unclosed "(" on line 4; the mistake is on line 3.
    with pytest.raises(ScriptParseError) as excinfo:
        parse_script("x = 1\ny = 2\nf = featurize(\nrank_on f\n")

    assert excinfo.value.line == 3
    assert excinfo.value.message == "'(' was never closed"


def test_positional_argument_error_downstream_of_a_call_that_closes_later_is_not_mistaken_for_unclosed():
    # The "(" on line 1 does close on line 4, so Lark's diagnosis of "y=," stands.
    with pytest.raises(ScriptParseError) as excinfo:
        parse_script("f = featurize(\n    x=1,\n    y=,\n)\nrank_on f\n")

    assert excinfo.value.line == 3
    assert excinfo.value.column == 7
    assert "never closed" not in excinfo.value.message
    assert "Unexpected token" in excinfo.value.message
