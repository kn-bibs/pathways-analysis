from patapy import render_text_table
from methods import Method, MethodResult


class MyMethod(Method):
    name = 'my analysis method'
    help = ''

    def run(self, experiment):
        pass


class MyResult(MethodResult):

    columns = ['score', 'p_value']


class ResultRow:

    def __init__(self, score, p_value):
        self.score = score
        self.p_value = p_value


def test_render_text_table(capsys):
    method = MyMethod()
    capsys.readouterr()     # clean buffers

    results = MyResult(
        [ResultRow(1, 0.1), ResultRow(2, 0.04)],
        files=['more-info.html'],
        description='Some important message'
    )

    render_text_table(method, results)
    std, err = capsys.readouterr()

    # message
    assert 'Some important message' in std

    # table
    for row in ['score	p_value', '1	0.1', '2	0.04']:
        assert row in std

    # files
    assert 'more-info.html' in std

    # no errors
    assert not err
