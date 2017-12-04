from contextlib import contextmanager
from tempfile import TemporaryFile

import pytest

from models import Gene, Sample, SampleCollection


@contextmanager
def temp_text_file(content):
    with TemporaryFile(mode='w+t') as file:
        file.write(content)
        file.seek(0)
        yield file


def test_init():
    genes1 = {Gene('BAD'): 1.2345, Gene('FUCA2'): 6.5432}
    genes2 = {Gene('BAD'): 2.3456, Gene('FUCA2'): 7.6543}

    samples = [Sample('Tumour_1', genes1), Sample('Tumour_2', genes2)]

    sample_collection = SampleCollection('Tumour', samples)

    assert sample_collection.name == 'Tumour'
    assert all(isinstance(k, Sample) for k in sample_collection.samples)


def test_from_csv():
    pass


gct_contents = """\
#1.2
2	3
NAME	DESCRIPTION	NORM-1	GBM-1	GBM-2	OV-1
TP53	na	348.61	172.52	236.45	130.2
MDM2	na	42.11	55.5	44.81	39.32
"""


def test_from_gct():
    with temp_text_file(gct_contents) as gct_file:

        with pytest.warns(
                UserWarning,
                match='Samples count (4) does not match with the 3 declared in All samples file.'
        ):
            collection = SampleCollection.from_gct_file('All samples', gct_file)

        assert len(collection.samples) == 4

        assert collection.labels == ['NORM-1', 'GBM-1', 'GBM-2', 'OV-1']
