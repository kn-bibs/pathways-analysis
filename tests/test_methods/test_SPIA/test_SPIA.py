from test_command_line.utilities import parse
from test_command_line.utilities import parsing_output


def test_help(capsys):
    with parsing_output(capsys) as text:
        parse('SPIA --help')
    assert 'signaling pathway impact analysis (SPIA)' in text.std


