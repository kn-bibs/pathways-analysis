from argparse import Namespace
from contextlib import contextmanager

import pytest

from command_line.main import SingleFileExperimentFactory
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


def test_slice():
    pass


def test_indices():
    pass


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
        't_2.tsv': (
            'Gene	Tumour_3	Tumour_4',
            'TP53	7	6',
            'BRCA2	6	8',
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

    def p_parse(command_line):
        """Parse with prefix"""
        return parse(prefix + command_line)

    accepted_commands = [
        'case t.tsv control c.tsv',
        'case t.tsv --samples Tumour_1 control c.tsv --samples Control_1',
        'case t.tsv --samples Tumour_1 control c.tsv --samples Control_1,Control_2',
        'data merged.tsv --case 1 --control 2',
        'data merged.tsv --case 2: --control :2',
        'control c.tsv case t.tsv t_2.tsv',
        # take all columns from t.tsv file and first column from t_2.tsv
        'control c.tsv case t.tsv t_2.tsv --columns 1,2 1'
    ]

    # TODO: this is not even a blackbox but will be developed soon
    for command in accepted_commands:
        print(command)
        p_parse(command)

    with pytest.raises(ValueError, message='columns for 2 files provided, expected for 1'):
        # the user should use --columns 1,2 instead
        p_parse('control c.tsv case t.tsv --columns 1 2')

    with pytest.raises(ValueError, message='columns for 1 files provided, expected for 2'):
        p_parse('control c.tsv case t.tsv t_2.tsv --columns 1')

    """
    cases = {
        # command => expected Namespace resulting from given command
        # all the other columns 5 >= are controls
        'case t.tsv --samples :4 --control c.tsv': '',
        'data merged.tsv --case :4': '',
        'data merged.tsv --control :3': '',
        'data merged.tsv --control 1,3' # all other are cases: '',
        
        # we want to ignore column number three altogether
        'data merged.tsv --case 1,2,4 --control 5-8': '',
        'data merged.tsv --case 1,2,4 --control 5:8': '',
        'control controls.ext --columns 1:4 case case.ext --case 1:3': '',
        
        # TODO: do we usually have sample names in headers?
        # dnusnf should only affect 'case':
        (
            'control controls.ext --columns 1:4'
            ' case case.ext --columns 1:3'
            ' --do_not_use_sample_names_from_header'
        ): '',
        (
            'control controls.ext --name "20_degrees" --columns 1:4'
            ' case case.ext --name "40_degrees" --columns 1:3'
        ): '',
        #'control controls.ext 1:4 case case.ext 1-3': '',
        # test multiple cases
        #'control controls.ext 1:4 case "20_degrees" case.ext 1-3 case "40_degrees" case.ext 4-5': '',
        #
        #'control controls.ext 1:4 case "20_degrees" case.ext 1-3 case "20_degrees" case2.ext 4-5': ''
    }
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


def test_general_help(capsys):

    with parsing_output(capsys) as text:
        parse('--help')

    for method in Method.members:
        assert method in text.std


def test_sub_parsers_help(capsys):
    # is custom sub-parser screen displayed and description used included in it?
    SingleFileExperimentFactory.description = 'A dummy description'
    SingleFileExperimentFactory.epilog = 'A dummy epilog'

    with parsing_output(capsys) as text:
        parse('data --help')

    lines = text.std.split('\n')
    half_of_output = len(lines) // 2

    # is description in the text and is it in the first 50% of lines?
    assert any(
        SingleFileExperimentFactory.description in line
        for line in lines[:half_of_output]
    )
    # is the epilog in the text and is it in the last 50% of lines?
    assert any(
        SingleFileExperimentFactory.epilog in line
        for line in lines[half_of_output:]
    )


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
