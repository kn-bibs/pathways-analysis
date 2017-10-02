import pytest
from test_command_line.utilities import parsing_output
from test_command_line.utilities import parse

from command_line.main import SingleFileExperimentFactory
from methods import Method
from models import Sample


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
            'TP53	6	6	7	6',
            'BRCA2	6	7	7	9',
        )
    }

    tmpdir.chdir()

    for filename, lines in files.items():
        file_object = tmpdir.join(filename)
        file_object.write('\n'.join(lines))


def p_parse(command_line):
    """Parse with prefix"""
    prefix = 'gsea '
    # as some method is always obligatory, each parser execution
    # will be prefixed with this method selection command
    return parse(prefix + command_line)


def test_files_loading(test_files):

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
        try:
            p_parse(command)
        except (Exception, SystemExit) as e:
            # we need the command to find out what caused the error
            print(command)
            raise e

    # TODO: make validation errors raise in parser context (but still easily testable!)

    with pytest.raises(ValueError, message='columns for 2 files provided, expected for 1'):
        # the user should use --columns 1,2 instead
        p_parse('control c.tsv case t.tsv --columns 1 2')

    with pytest.raises(ValueError, message='columns for 1 files provided, expected for 2'):
        p_parse('control c.tsv case t.tsv t_2.tsv --columns 1')

    with pytest.raises(ValueError, message='Cannot handle data and case/control at once'):
        p_parse('data merged.tsv control c.tsv')

    with pytest.raises(ValueError, message='Neither data nor (case & control) have been provided!'):
        p_parse('')


def test_columns_purpose_deduction(test_files):

    def make_samples(samples_dict):
        return [
            Sample.from_names(name, values)
            for name, values in samples_dict.items()
        ]

    # all the other columns (id >= 2) are cases
    commands = [
        'data merged.tsv --control :2',
        'data merged.tsv --control 0,1'
    ]
    expected_cases = {
        'Tumour_1': {'TP53': 7, 'BRCA2': 7},
        'Tumour_2': {'TP53': 6, 'BRCA2': 9},
    }
    expected_controls = {
        'Control_1': {'TP53': 6, 'BRCA2': 6},
        'Control_2': {'TP53': 6, 'BRCA2': 7},
    }
    for command in commands:

        print(command)
        opts = p_parse(command)

        assert opts.control.phenotype.samples == make_samples(expected_controls)
        assert opts.case.phenotype.samples == make_samples(expected_cases)


def TODO():
    """
    cases = {
        # command => expected Namespace resulting from given command
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

    # if there are no arguments provided, the parser should
    # show the usage summary (and do not raise any errors)
    with parsing_output(capsys) as text:
        parse(None)

    assert 'usage' in text.err


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
