import re

from multiprocess.progress_bar import progress_bar


def test_bar(capfd):
    capfd.readouterr()     # clean buffers

    data = [1, 1, 1, 1, 1]

    with progress_bar(len(data)) as queue:
        # simulate five consecutive updates
        for i in data:
            queue.put(i)

    assert queue.empty()

    # let's capture the output on both standard stream
    std, err = capfd.readouterr()

    # at the end there should be 5 iterations observed, the timing does not matter
    expected = re.compile(
        r'(.*?)100%|##########| 5/5 \[\d\d:\d\d<\d\d:\d\d, \d+\.\d\dit/s\]\n$',
        re.MULTILINE
    )

    # tqdm writes to err stream
    assert expected.match(err)
