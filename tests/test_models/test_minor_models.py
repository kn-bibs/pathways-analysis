from models import Gene, Sample

# TODO method for creating samples faster


def test_gene_init():
    gene = Gene('BAD')
    assert gene.name == 'BAD'

    same = Gene('BAD')
    assert same == gene
    assert same is gene

    g = Gene('TP53', 'Tumour suppressor p53')
    assert g.name == 'TP53'
    assert g.description == 'Tumour suppressor p53'

    from copy import copy
    assert copy(g).id is g.id


def test_sample_init():
    genes = {Gene('BAD'): 1.2345, Gene('FUCA2'): 6.5432}

    sample = Sample('Tumour_1', genes)

    assert sample.name == 'Tumour_1'
    assert all(isinstance(k, Gene) for k in sample.data.keys())
    assert sample.data == genes


def test_sample_from_names():
    data = {'BAD': 1.2345, 'FUCA2': 6.5432}

    sample = Sample.from_names('Tumour_1', data)

    assert sample.name == 'Tumour_1'
    assert all(isinstance(k, Gene) for k in sample.data.keys())
    assert [k.name for k in sample.data.keys()] == ['BAD', 'FUCA2']
    assert all(isinstance(v, float) for v in sample.data.values())
    assert list(sample.data.values()) == [1.2345, 6.5432]

