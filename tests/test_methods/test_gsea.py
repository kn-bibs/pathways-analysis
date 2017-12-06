from test_command_line.utilities import parse
from test_command_line.utilities import parsing_output


def test_database(capsys):
    with parsing_output(capsys) as text:
        parse('gsea database --show_gene_sets')
    assert 'HALLMARK_HYPOXIA' in text.std


def test_licence(capsys):
    with parsing_output(capsys) as text:
        parse('gsea --show_licence')
    assert 'copyright (c) 2004-2017 Broad Institute, Inc.' in text.std
