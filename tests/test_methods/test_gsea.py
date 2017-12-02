from test_command_line.utilities import parse
from test_command_line.utilities import parsing_output

from methods.gsea import GeneralisedGSEA
from methods.gsea.gsea import ScoreDistribution
from methods.gsea.metrics import DifferenceOfClasses
from methods.gsea.signatures import MolecularSignatureDatabase, GeneSet
from models import SampleCollection, Sample, Gene


def create_test_db():
    sets = [
        GeneSet('anthrax pathway', ['MAP2K2', 'MAP2K1']),
        GeneSet('p53 pathway', ['MDM2', 'TP53'])
    ]
    return MolecularSignatureDatabase({
        s.name: s for s in sets
    })


def test_database(capsys):
    with parsing_output(capsys) as text:
        parse('gsea database --show_gene_sets')
    assert 'HALLMARK_HYPOXIA' in text.std

    db = create_test_db()
    assert len(db.gene_sets) == 2


def test_gene_set():
    gene_set = GeneSet('set', ['MDM2', 'TP53'])
    assert Gene('TP53') in gene_set


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
    assert GeneralisedGSEA.estimate_significance_level(9, ScoreDistribution([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])) == 0.1
    # TODO: possibly this could return LessThan
    assert GeneralisedGSEA.estimate_significance_level(9, ScoreDistribution([-1, -2, -3])) == 0


def test_ranked_list():
    db = create_test_db()

    tp53 = Gene('TP53')
    map2k1 = Gene('MAP2K1')

    case = SampleCollection(
        'case',
        [Sample('1', {tp53: 2, map2k1: 1})]
    )
    control = SampleCollection(
        'control',
        [Sample('1', {tp53: 1, map2k1: 1})]
    )

    gsea = GeneralisedGSEA(database=db, ranking_metric=DifferenceOfClasses)

    assert gsea.create_ranked_gene_list(case, control) == [(tp53, 1), (map2k1, 0)]
