import pytest

from models import Experiment, SampleCollection, Sample

# TODO methods for creating samples and sample_collections


def test_experiment_init():
    data1 = {'BAD': 1.2345, 'FUCA2': 6.5432}
    data2 = {'BAD': 2.3456, 'FUCA2': 7.6543}
    data3 = {'BAD': 3.4567}

    tumour_samples = [Sample.from_names('Tumour_1', data1), Sample.from_names('Tumour_2',data2)]
    normal_samples = [Sample.from_names('Normal_1',data3)]

    tumour = SampleCollection('Tumour', tumour_samples)
    normal = SampleCollection('Normal', normal_samples)

    experiment = Experiment(case=tumour, control = normal)

    assert isinstance(experiment.case, SampleCollection)
    assert isinstance(experiment.control, SampleCollection)

    assert experiment.case == tumour
    assert experiment.control == normal

def test_experiment_from_tsv():
    pass

def test_experiment_from_gsea_file():
    pass

def test_experiment_get_all():
    pass

def test_experiment_fold_change():
    pass


