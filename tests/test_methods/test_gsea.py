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


def test_help(capsys):
    with parsing_output(capsys) as text:
        parse('gsea --help')
    assert 'GSEA method' in text.std


def test_more_extreme():
    from methods.gsea.gsea import is_more_extreme
    for a, b in [(10, 5), (-10, -5)]:
        assert is_more_extreme(a, b)
        assert not is_more_extreme(b, a)


def test_estimate_significance():
    from methods.gsea import GeneralisedGSEA
    from methods.gsea.gsea import ScoreDistribution
    assert GeneralisedGSEA.estimate_significance_level(9, ScoreDistribution([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])) == 0.1
    # TODO: this could return LessThan
    assert GeneralisedGSEA.estimate_significance_level(9, ScoreDistribution([-1, -2, -3])) == 0
