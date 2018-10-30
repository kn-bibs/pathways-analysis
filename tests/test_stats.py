from models import Sample, SampleCollection, Experiment
import pandas as pd
from stats import ttest


def test_ttest():
    data1 = {'BAD': 1.2345, 'FUCA2': 6.5432}
    data2 = {'BAD': 2.3456, 'FUCA2': 7.6543}
    data3 = {'BAD': 6.3456, 'FUCA2': 11.6543}
    data4 = {'BAD': 7.1111, 'FUCA2': 9.9711}

    tumour_samples = [Sample.from_names('Tumour_1', data1), Sample.from_names('Tumour_2', data2)]
    normal_samples = [Sample.from_names('Normal_1', data3), Sample.from_names('Normal_2', data4)]

    tumour = SampleCollection('Tumour', tumour_samples)
    normal = SampleCollection('Normal', normal_samples)

    experiment = Experiment(case=tumour, control=normal)
    tt = ttest(experiment)
    assert isinstance(tt, pd.Series)
    assert all(gene in list(tt.keys()) for gene in experiment.get_all().genes)
