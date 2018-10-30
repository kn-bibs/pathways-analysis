from databases import KEGGPathways
import networkx as nx
import pytest


def test_keggpathways_init():
    db = KEGGPathways()
    assert db.organism == "hsa"
    db1 = KEGGPathways("Gallus gallus")
    db2 = KEGGPathways("gallus gallus")
    db3 = KEGGPathways("chicken")
    db4 = KEGGPathways("Chicken")
    db5 = KEGGPathways("ChIcKeN")
    assert all([d.organism == "gga" for d in [db1, db2, db3, db4, db5]])


def test_search_by_gene():
    db = KEGGPathways()
    assert isinstance(db.search_by_gene('BRCA2'), dict)
    pathways = db.search_by_gene('TheMostImportantGene')
    assert pathways == {}


def test_get_pathway():
    db = KEGGPathways()
    pathway = db.get_pathway('hsa04630')
    assert all([pathway.node[node]['type'] == 'gene' for node in pathway.nodes])
    assert all(isinstance(attr, list) for attr in nx.get_edge_attributes(pathway, 'type').values())
    assert all([pathway.degree(node) > 0 for node in pathway.nodes])


def test_get_organism_code():
    with pytest.raises(KeyError):
        KEGGPathways().get_organism_code('Homo bioinformaticus')


def test_get_gene_code():
    gen = KEGGPathways().get_gene_code('BIOINF2018')
    assert gen == ''
