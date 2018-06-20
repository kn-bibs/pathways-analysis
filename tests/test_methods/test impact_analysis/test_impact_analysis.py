import networkx as nx
from methods.impact_analysis import *
from models import *


def minimal_data():
    brca2 = Gene('BRCA2')
    kit = Gene('KIT')

    case = SampleCollection(
        'case',
        [Sample('1', {brca2: 2, kit: 1})]
    )
    control = SampleCollection(
        'control',
        [Sample('2', {brca2: 1, kit: 1})]
    )
    experiment = Experiment(case, control)
    return experiment


def minimal_pathway():
    pathway = nx.DiGraph()
    genes = ['BRCA2', 'KIT', 'BRCA1']
    for i in range(len(genes)):
        pathway.add_node(genes[i], type=["gene"])
    edges = {('BRCA2', 'BRCA1'): 'activation',
             ('BRCA1', 'KIT'): 'repression'}
    for e in edges.keys():
        pathway.add_edge(e[0], e[1], type=[edges[e]])
    return pathway


def test_caluculate_corrections():
    pvalues = pd.Series([0.2, 0.1, 1])
    ia = ImpactAnalysis()
    fdr, bonf = ia.calculate_corrections(pvalues)
    assert all([x is not None for x in [fdr, bonf]])
    assert all([x.size == len(pvalues) for x in [fdr, bonf]])


def test_calculate_impact_factor():
    ia = ImpactAnalysis()
    experiment = minimal_data()
    ia.FC = experiment.calculate_fold_change()
    ia.degs = pd.Series({Gene('BRCA2'): True})
    ia.experiment_genes = set([x.name for x in experiment.get_all().genes])
    pathway = minimal_pathway()
    pathway_if, pvalue = ia.calculate_impact_factor(experiment, pathway)

    assert pathway_if == 7.5
    assert pvalue == 1.0

    pathway.add_edge("KIT", "BRCA2", type=["activation"])

    pathway_if, pvalue = ia.calculate_impact_factor(experiment, pathway)

    assert pathway_if == 10.5
    assert pvalue == 1.0
