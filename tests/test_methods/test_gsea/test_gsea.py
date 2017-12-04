import random

import numpy
from test_command_line.utilities import parse
from test_command_line.utilities import parsing_output

from methods.gsea import GeneralisedGSEA
from methods.gsea.gsea import ScoreDistribution
from methods.gsea.metrics import difference_of_classes
from methods.gsea.signatures import MolecularSignatureDatabase, GeneSet
from models import SampleCollection, Sample, Gene, Experiment


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
    assert str(gene_set) == '<GeneSet: set with 2 genes>'


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


def test_score_distribution():
    dist = ScoreDistribution([-1, 1])
    assert dist.negative_scores == [-1]
    assert dist.positive_scores == [+1]
    assert len(dist) == 2


def test_estimate_significance():
    assert GeneralisedGSEA.estimate_significance_level(
        9,
        ScoreDistribution([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    ) == 0.1
    # TODO: possibly this could return LessThan
    assert GeneralisedGSEA.estimate_significance_level(
        9,
        ScoreDistribution([-1, -2, -3])
    ) == 0


def minimal_data():

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

    return tp53, map2k1, case, control


def test_ranked_list():
    db = create_test_db()
    tp53, map2k1, case, control = minimal_data()

    gsea = GeneralisedGSEA(database=db, ranking_metric=difference_of_classes)

    assert gsea.create_ranked_gene_list(case, control) == [(tp53, 1), (map2k1, 0)]


def test_shuffle_samples():
    from methods.gsea.shufflers import shuffle_and_divide

    tp53, map2k1, case, control = minimal_data()

    merged_collection = case + control

    numpy.random.seed(0)
    new_case, new_control = shuffle_and_divide(merged_collection, 1)

    assert new_case.samples == control.samples
    assert new_control.samples == case.samples


def test_run():
    tp53, map2k1, case, control = minimal_data()
    experiment = Experiment(case, control)

    gsea = GeneralisedGSEA(create_test_db(), min_genes=3)
    gene_sets = gsea.trim_gene_sets(gsea.gene_sets, experiment)
    assert not gene_sets

    db = create_test_db()

    # normalized enrichment score is depends null distribution,
    # which is a result of random shuffling; setting seed to
    # have reproducible results
    numpy.random.seed(0)
    random.seed(0)

    gsea = GeneralisedGSEA(
        db,
        ranking_metric=difference_of_classes,
        min_genes=1,
        processes=1
    )

    results = gsea.run(experiment)
    assert len(gsea.gene_sets) == 2

    assert results.scored_list == [
        db.gene_sets['anthrax pathway'],
        db.gene_sets['p53 pathway']
    ]

    anthrax = results.scored_list[0]

    # the enrichment is set to zero to reflect the fact of no
    # change in expression with respect to given metric (rank
    # for difference of classes: 1 - 1 = 0)
    assert anthrax.enrichment == 0
    assert anthrax.nominal_p_value == 0

    p53 = results.scored_list[1]

    # caveat: this is hardened (not hand-calculated) result;
    # would be beneficial to try to hand-calculate this too.
    assert p53.enrichment == 1.4634615384615386
