from contextlib import contextmanager
from tempfile import TemporaryFile

from pytest import warns

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


csv_contents = """\
NAME,NORM-1,GBM-1,GBM-2,OV-1
TP53,348.61,172.52,236.45,130.2
MDM2,42.11,55.5,44.81,39.32
"""


def test_from_csv():

    with temp_text_file(csv_contents) as csv_file:

        collection = SampleCollection.from_csv_file('all_samples.csv', csv_file)

        assert len(collection.samples) == 4

    with temp_text_file(csv_contents) as csv_file:

        with warns(UserWarning, match='You are using not comma delimiter for what looks like csv file.'):
            SampleCollection.from_csv_file('all_samples.csv', csv_file, delimiter='\t')


# note: number of samples defined incorrectly (3 instead of four) on purpose
gct_contents = """\
#1.2
2	3
NAME	DESCRIPTION	NORM-1	GBM-1	GBM-2	OV-1
TP53	na	348.61	172.52	236.45	130.2
MDM2	na	42.11	55.5	44.81	39.32
"""

# it seems that the GCT files from the Broad Institute's
# Cancer Cell Line Encyclopedia (CCLE) have an additional
# tab character at the end of the version and counts lines;
# this additional tab needs to be handled properly.
ccle_like_gct = """\
#1.2	
1	4	
NAME	DESCRIPTION	NORM-1	GBM-1	GBM-2	OV-1
TP53	na	348.61	172.52	236.45	130.2
"""


def test_from_gct():

    with temp_text_file(gct_contents) as gct_file:

        expected_warning = 'Samples count \(4\) does not match with the 3 declared in all_samples.gct file.'

        with warns(UserWarning, match=expected_warning):
            collection = SampleCollection.from_gct_file('all_samples.gct', gct_file)

        assert len(collection.samples) == 4

        assert collection.labels == ['NORM-1', 'GBM-1', 'GBM-2', 'OV-1']

    # replace version definition
    lines = gct_contents.split('\n')
    lines[0] = '#1.1'
    old_content = '\n'.join(lines)

    with temp_text_file(old_content) as old_gct_file:
        with warns(UserWarning, match='Unsupported version of GCT file'):
            SampleCollection.from_gct_file('Outdated file', old_gct_file)

    with temp_text_file(ccle_like_gct) as gct_file:
        with warns(None) as warnings:
            collection = SampleCollection.from_gct_file('all_samples.gct', gct_file)
        assert not warnings.list
