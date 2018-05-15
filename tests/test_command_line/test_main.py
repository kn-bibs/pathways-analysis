import pytest
from pytest import fixture

from command_line.main import SingleFileExperimentFactory, CLI
from methods import Method
from models import Sample, Gene

from .utilities import parsing_error
from .utilities import parsing_output
from .utilities import parse


def make_samples(samples_dict):
    """Create samples from dict representation"""
    return [
        Sample.from_names(name, values)
        for name, values in samples_dict.items()
    ]


expected_cases = make_samples({
    'Tumour_1': {'TP53': 7, 'BRCA2': 7},
    'Tumour_2': {'TP53': 6, 'BRCA2': 9},
})

expected_controls = make_samples({
    'Control_1': {'TP53': 6, 'BRCA2': 6},
    'Control_2': {'TP53': 6, 'BRCA2': 7},
})


@fixture
def test_files(tmpdir):
    # Here I assume that all the files are in TSV format (for now, we may want to change this in the future)

    # create temporary files
    files = {
        # the ".ext" extensions are here just to visually
        # mark that the strings represent some file names.
        'c.tsv': (
            'Gene	Control_1	Control_2',
            'TP53	6	6',
            'BRCA2	6	7',
        ),
        't.tsv': (
            'Gene	Tumour_1	Tumour_2',
            'TP53	7	6',
            'BRCA2	7	9',
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
    create_files(tmpdir, files)


def create_files(tmpdir, files):

    tmpdir.chdir()

    for filename, lines in files.items():
        file_object = tmpdir.join(filename)
        file_object.write('\n'.join(lines))


class DummyMethod(Method):
    name = 'dummy'
    help = ''

    def run(self, experiment):
        pass


def p_parse(command_line):
    """Parse with prefix"""
    prefix = 'dummy '
    # as some method is always obligatory, each parser execution
    # will be prefixed with this method selection command
    try:
        return parse(prefix + command_line)
    except (Exception, SystemExit) as e:
        # we need the command to find out what caused the error
        print(command_line)
        raise e


def test_simple_files_loading(test_files):
    """Only general tests, no details here"""

    # one file for case, one for control
    opts = p_parse('case t.tsv control c.tsv')
    assert len(opts.case.sample_collection.samples) == 2

    # two files for case, one for control
    opts = p_parse('control c.tsv case t.tsv t_2.tsv')
    assert len(opts.case.sample_collection.samples) == 4

    with parsing_error(match='Neither data nor \(case & control\) have been provided!'):
        p_parse('')

    sample_collections = {'case': 'Control', 'control': 'Case'}
    for sample_collection, name in sample_collections.items():
        with parsing_error(match=f'{name} has not been provided!'):
            p_parse(f'{sample_collection} c.tsv')


def test_select_samples(test_files):

    # lets select only first samples from both files
    commands = [
        'case t.tsv --columns 0 control c.tsv --columns 0',
        'case t.tsv --samples Tumour_1 control c.tsv --samples Control_1'
    ]
    for command in commands:
        opts = p_parse(command)

        assert opts.control.sample_collection.samples == expected_controls[:1]
        assert opts.case.sample_collection.samples == expected_cases[:1]

    # get both tumour samples from file t.tsv and
    # the first sample (Tumour_3) from file t_2.tsv
    commands = [
        'control c.tsv case t.tsv t_2.tsv --columns 0,1 0',
        'case t.tsv t_2.tsv --samples Tumour_1,Tumour_2 Tumour_3 control c.tsv'
    ]

    for command in commands:
        opts = p_parse(command)
        assert len(opts.case.sample_collection.samples) == 3

    # lets try to grab a sample which is not in the file

    expected_message = (
        "Samples {'Control_1'} are not available in t.tsv file.\n"
        "Following samples were found: Tumour_1, Tumour_2."
    )

    with parsing_error(match=expected_message):
        p_parse('case t.tsv --samples Control_1 control c.tsv')

    with parsing_error(match='columns for 2 files provided, expected for 1'):
        # the user should use --columns 0,1 instead
        p_parse('control c.tsv case t.tsv --columns 0 1')

    with parsing_error(match='columns for 1 files provided, expected for 2'):
        p_parse('control c.tsv case t.tsv t_2.tsv --columns 1')


def test_columns_purpose_deduction(test_files):

    # all the other columns (id >= 2) are cases
    commands = [
        'data merged.tsv --control :2',
        'data merged.tsv --control 0,1',
        'data merged.tsv --case 2:',
        'data merged.tsv --case 2,3'
    ]
    for command in commands:
        opts = p_parse(command)

        assert opts.control.sample_collection.samples == expected_controls
        assert opts.case.sample_collection.samples == expected_cases


def test_non_tab_delimiter(tmpdir):

    create_files(tmpdir, {
        'c.tsv': (
            'Gene,Control_1,Control_2',
            'TP53,6,6',
            'BRCA2,6,7',
        ),
        't.tsv': (
            'Gene	Tumour_1	Tumour_2',
            'TP53	7	6',
            'BRCA2	7	9',
        ),
    })

    opts = p_parse('case t.tsv control c.tsv --delimiter ,')

    assert opts.control.sample_collection.samples == expected_controls
    assert opts.case.sample_collection.samples == expected_cases


def test_file_with_description(test_files, tmpdir):

    create_files(tmpdir, {
        'control_with_descriptions.tsv': (
            'Gene	Description	Control_1	Control_2',
            'TP53	Tumour protein 53	6	6',
            'BRCA2	Breast cancer type 2 s. protein	6	7',
        )
    })

    expected_warning = (
        'First line of your file contains "description" column, '
        'but you did not provide "--description_column" argument.'
    )

    # user forgot
    with pytest.warns(UserWarning, match=expected_warning):
        opts = p_parse('case t.tsv control control_with_descriptions.tsv')
        assert len(opts.control.sample_collection.samples) == 3

    # user remembered
    opts = p_parse('case t.tsv control control_with_descriptions.tsv -d')
    assert len(opts.control.sample_collection.samples) == 2

    assert set(opts.control.sample_collection.samples[0].genes) == {
        Gene('TP53', description='Tumour protein 53'),
        Gene('BRCA2', description='Breast cancer type 2 s. protein')
    }


def test_custom_sample_names(test_files):

    opts = p_parse(
        'case t.tsv control c.tsv t_2.tsv --header my_control their_control'
    )

    # case should not be affected anyhow there
    assert opts.case.sample_collection.samples == expected_cases

    controls = opts.control.sample_collection.samples

    # are two files loaded? (each have two samples)
    assert len(controls) == 4

    sample_names = {control.name for control in controls}

    expected_names = {
        'my_control_1', 'my_control_2',
        'their_control_1', 'their_control_2'
    }

    assert expected_names == sample_names


def test_merged_file(test_files):
    # advanced columns purpose inferring/deduction is tested separately
    # in `test_columns_purpose_deduction`

    # merged.tsv has two controls and two tumours, in this order

    commands_first_samples = [
        'data merged.tsv --case 2 --control 0',
        'data merged.tsv --case 2,2 --control 0,0',
    ]
    for command in commands_first_samples:
        opts = p_parse(command)

        assert opts.control.sample_collection.samples == expected_controls[:1]
        assert opts.case.sample_collection.samples == expected_cases[:1]

    commands_all_samples = [
        'data merged.tsv --case 2: --control :2',
        'data merged.tsv --case 2,3 --control 0:2',
        'data merged.tsv --case 2-4 --control 0-2'
    ]
    for command in commands_all_samples:
        opts = p_parse(command)

        assert opts.control.sample_collection.samples == expected_controls
        assert opts.case.sample_collection.samples == expected_cases

    with parsing_error(match='Neither --case nor --control provided'):
        p_parse('data merged.tsv')

    with parsing_error(match='Cannot handle data and case/control at once'):
        p_parse('data merged.tsv --case 1 --control 2 control c.tsv')


def test_general_help(capsys):

    with parsing_output(capsys) as text:
        parse('--help')

    for method in Method.members:
        assert method in text.std

    with parsing_error(match='unrecognized arguments: --controll'):
        p_parse('data merged.tsv --case 1 --controll 2 control c.tsv')


def test_shows_usage_when_no_args(capsys):

    # if there are no arguments provided, the parser should
    # show the usage summary (and do not raise any errors)
    with parsing_output(capsys) as text:
        parse(None)

    assert 'usage' in text.err


def test_sub_parsers_help(capsys):
    # do we get the `name` substitution right?
    SingleFileExperimentFactory.__doc__ = 'Description of {parser_name}'
    cli = CLI()
    assert cli.all_subparsers['data'].parser_name == 'data'

    with parsing_output(capsys) as text:
        parse('data --help')

    assert 'Description of data' in text.std

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
