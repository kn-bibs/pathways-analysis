import pytest

from methods import Method
from command_line import CLI


# Here I assume that all the files are in TSV format (for now, we may want to change this in the future)


def parse(command_line):
    """Helper executing Argument parser

    Returns:
        Namespace with parsed arguments
    """
    commands = command_line.split(' ')
    return CLI().parse(commands)


def test_files_loading(tmpdir):
    # method is obligatory, each parser execution will
    # be prefixed with this method selection command

    prefix = '--method gsea '

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

    p = {}

    for filename, lines in files.items():
        file_object = tmpdir.join(filename)
        p[filename] = str(file_object)
        file_object.write('\n'.join(lines))

    accepted_commands = [
        f'case --files {p["t.tsv"]} control --files {p["c.tsv"]}',
        f'case --files {p["t.tsv"]} --samples Tumour_1 control --files {p["c.tsv"]} --samples Control_1',
        f'case --files {p["t.tsv"]} --samples Tumour_1 control --files {p["c.tsv"]} --samples Control_1 Control_2',
        #f'data --files {p["merged.tsv"]} --case 1 --control 2'
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


def test_methods(capsys):
    # TODO: needs implementation
    help_command = '--methods -h'
    methods_help = parse(help_command)

    for method in Method.members:
        assert method.name in methods_help

    class MyMethod(Method):

        name = 'some_method'
        help = 'Some important text help'

        def __init__(self, my_argument=None):
            super().__init__()
            if my_argument == 'value':
                raise Success('Correct Argument provided +1')

        def run(self, experiment):
            raise Success('Run +10')

    # test Argument
    with pytest.raises(Success):
        parse('--method some_method --my_argument value')

    # test Argument
    with pytest.raises(Success):
        parse('--method some_method --my_argument value')

    # test help
    parse('--method some_method --help')
    out, _ = capsys.readouterr()

    assert 'Some important text help' in out
