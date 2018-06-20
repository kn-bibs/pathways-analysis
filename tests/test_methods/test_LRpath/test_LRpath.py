from test_command_line.utilities import parse
from test_command_line.utilities import parsing_output

from methods.LRpath import LRpath

from models import Experiment, Gene, SampleCollection, Sample


def test_help(capsys):
    with parsing_output(capsys) as text:
        parse('LRpath --help')
    assert 'LRpath performs gene' in text.std


def create_test_db():
    db = {'GO:0001101': ['8649'], 'GO:0001819': ['6288', '3600']}
    return db


def minimal_data():

    tp53 = Gene('TP53')
    map2k1 = Gene('MAP2K1')

    case = SampleCollection(
        'case',
        [Sample('1', {tp53: 4, map2k1: 5})]
    )
    control = SampleCollection(
        'control',
        [Sample('1', {tp53: 1, map2k1: 1})]
    )

    return tp53, map2k1, case, control


def test_run():

    tp53, map2k1, case, control = minimal_data()
    experiment = Experiment(case, control)
    db = create_test_db()
    lrpath = LRpath(min_g=1, database=db)

    results = lrpath.run(experiment=experiment)
    assert len(results.scored_list) == 2

    match = experiment.case.as_array()
    data, names = lrpath.create_data(match)
    data2, geneid = lrpath.name_geneid(data, experiment.case.genes)
    assert len(geneid) == 2
    assert geneid == ['hsa:6288', 'hsa:8649']

    assert len(data2) == 2
    assert data2.index.all(geneid)

    tp53 = results.scored_list[0]
    assert round(tp53.LRcoeff, 10) == round(0.03968345542274243, 10)
    assert round(tp53.LRpvalue, 10) == round(0.9124653704749152, 10)
    assert round(tp53.odds_ratio, 10) == round(0.7814398312529739, 10)

    map2k1 = results.scored_list[1]
    assert round(map2k1.LRcoeff, 10) == round(-0.0396834554227424, 10)
