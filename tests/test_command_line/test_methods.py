import pytest

from declarative_parser.parser import Argument
from patapy import run
from methods import Method
from test_command_line.test_main import test_files
from test_command_line.utilities import parse
from test_command_line.utilities import parsing_output


class Success(Exception):
    pass


def test_methods(capsys, test_files):
    suffix = 'case t.tsv control c.tsv'

    class MyMethod(Method):

        name = 'some_method'
        help = 'Some important text help'

        other_argument = Argument(
            type=int,
            help='Only integer numbers!',
            default=0,
        )

        def __init__(self, mode, my_argument: float=None, other_argument=None):
            """

            Args:
                mode: This argument has to be passed!
                my_argument: Help and documentation for my_argument
                other_argument: Documentation for the other_argument
            """

            if mode == 'active':
                raise Success('In active mode!')
            if my_argument == 1.4:
                raise Success('Correct Argument provided +1')
            if type(other_argument) is int:
                raise Success(f'Parsed +{other_argument}')

        def run(self, experiment):
            raise Success('Run +10')

    # test mode which is a required argument
    with pytest.raises(Success, message='In active mode!'):
        parse(f'some_method active {suffix}')

    # test my_argument
    with pytest.raises(Success, message='Correct Argument provided +1'):
        parse(f'some_method a --my_argument 1.4 {suffix}')

    # test run
    with pytest.raises(Success, message='Run +10'):
        run(f'run.py some_method a --my_argument 2.1 {suffix}'.split())

    # test other_argument
    with pytest.raises(Success, message='Parsed +5'):
        parse(f'some_method a --other_argument 5 {suffix}')

    # test help
    for test_command in ['some_method --help', '-h some_method']:
        with parsing_output(capsys, contains='Some important text help'):
            parse(test_command)

    with parsing_output(capsys, does_not_contain='Some important text help'):
        parse('--help')

    # test arguments help
    arguments_help = {
        'This argument has to be passed!': True,
        'Help and documentation for my_argument': True,
        'Only integer numbers!': True,
        # Following should not be shown, as it is overwritten by
        # "integers-only" text above from Argument definition:
        'Documentation for the other_argument': False
    }

    for text, is_expected in arguments_help.items():
        contains, does_not_contain = None, None
        if is_expected:
            contains = text
        else:
            does_not_contain = text

        with parsing_output(capsys, contains=contains, does_not_contain=does_not_contain):
            parse('some_method -h')

