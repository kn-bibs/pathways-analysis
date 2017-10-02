import pytest
from test_command_line.utilities import parsing_output
from test_command_line.utilities import parse

from command_line.main import SingleFileExperimentFactory
from methods import Method


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
