from argparse import Namespace
from contextlib import contextmanager

import pytest

from command_line.parser import Argument
from methods import Method
from command_line import CLI
from patapy import run


def parse(command_line):
    """Helper executing Argument parser

    Returns:
        Namespace with parsed arguments
    """
    commands = command_line.split(' ')
    return CLI().parse(commands)


@pytest.fixture
def test_files(tmpdir):
    # Here I assume that all the files are in TSV format (for now, we may want to change this in the future)

    # create temporary files
    files = {
        # the ".ext" extensions are here just to visually
        # mark that the strings represent some file names.
        'c.tsv': (
            'Gene	Control_1	Control_2',
            'TP53	1	2',
            'BRCA2	4	5',
        ),
        't.tsv': (
            'Gene	Tumour_1	Tumour_2',
            'TP53	6	6',
            'BRCA2	6	7',
        ),
        'control_without_headers.tsv': (
            'TP53	5	6',
            'BRCA2	4	8',
        ),
        'merged.tsv': (
            'Gene	Control_1	Control_2	Tumour_1	Tumour_2',
            'TP53	6	6	6	6',
            'BRCA2	6	7	6	7',
        )
    }

    tmpdir.chdir()

    for filename, lines in files.items():
        file_object = tmpdir.join(filename)
        file_object.write('\n'.join(lines))


def test_files_loading(test_files):
    # method is obligatory, each parser execution will
    # be prefixed with this method selection command

    prefix = 'gsea '

    accepted_commands = [
        'case --files t.tsv control --files c.tsv',
        'case --files t.tsv --samples Tumour_1 control --files c.tsv --samples Control_1',
        'case --files t.tsv --samples Tumour_1 control --files c.tsv --samples Control_1,Control_2',
        'data --files merged.tsv --case 1 --control 2'
    ]

    # TODO: this is not even a blackbox but will be developed soon
    for command in accepted_commands:
        print(command)
        parse(prefix + command)

    """
    cases = {
        # command => expected Namespace resulting from given command
        '': '',
        '--control control.ext --case tumour_1.ext': '',
        '--control some_file.ext --case 20_degrees.ext 40_degrees.ext': '',
        # all the other columns 5 >= are controls
        '--case some_file.ext --case_samples :4 --control som_file.ext': '',
        '--data some_file.ext --case_samples :4': '',
        '--data some_file.ext --control_samples :3': '',
        '--data some_file.ext --control_samples 1,3 # all other are cases': '',
        # we want to ignore column number three altogether
        '--data some_file.ext --case_samples 1,2,4 --control_samples 5-8': '',
        '--data some_file.ext --case_samples 1,2,4 --control_samples 5:8': '',
        '--control controls.ext --control_columns 1:4 --case case.ext --case_columns 1:3': '',
        # TODO: czy większość danych ma headery? Czy switch powinien być domyślnie włączony?'
        (
            '--control controls.ext --control_samples 1:4'
            ' --case case.ext --case_samples 1:3'
            ' --do_not_use_sample_names_from_header'
        ): '',
        (
            '--control controls.ext --name "20_degrees" --control_samples 1:4'
            ' --case case.ext --name "40_degrees" --case_samples 1:3'
        ): '',
        #' --control controls.ext 1:4 --case case.ext 1-3': '',
        # test multiple cases
        #'--control controls.ext 1:4 --case "20_degrees" case.ext 1-3 --case "40_degrees" case.ext 4-5': '',
        #
        #'--control controls.ext 1:4 --case "20_degrees" case.ext 1-3 --case "20_degrees" case2.ext 4-5': ''
    }
    for command, result in cases.items():
        assert parse(prefix + command) == result
    """


class Success(Exception):
    pass


@contextmanager
def parsing_output(capsys, contains=None, does_not_contain=None):
    text = Namespace()
    capsys.readouterr()     # clean buffers
    with pytest.raises(SystemExit):
        yield text
    text.std, text.err = capsys.readouterr()
    if contains:
        assert contains in text.std
    if does_not_contain:
        assert does_not_contain not in text.std


def test_help(capsys):

    with parsing_output(capsys) as text:
        parse('--help')

    for method in Method.members:
        assert method in text.std


def test_methods(capsys, test_files):
    suffix = 'case --files t.tsv control --files c.tsv'

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


def test_analyze_docstring():

    docstring = """Some docstring.
    
    Arguments:
        my_arg: is an important argument
        active: should some feature be active
                or maybe it should be not?
        spread:
            should be big or small?
            how big or how small?
            
    Example:
        examples should not be interpreted as an argument
    
    Returns:
        results
    """
    from command_line.method_parser import analyze_docstring
    args = analyze_docstring(docstring)

    expected_args = {
        'my_arg': 'is an important argument',
        'active': 'should some feature be active or maybe it should be not?',
        'spread': 'should be big or small? how big or how small?'
    }

    for name, value in expected_args.items():
        assert args[name] == value

    assert 'Example' not in args
