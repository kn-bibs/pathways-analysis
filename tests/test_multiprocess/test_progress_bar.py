import re

from multiprocess import STOP
from multiprocess.progress_bar import progress_bar, progress_bar_worker

# simulate five consecutive, step-by-step updates
data = [1, 1, 1, 1, 1]

# at the end there should be 5 iterations observed, the timing does not matter
expected_bar = re.compile(
    r'(.*?)100%|##########| 5/5 \[\d\d:\d\d<\d\d:\d\d, \d+\.\d\dit/s\]\n$',
    re.MULTILINE
)


def test_bar_worker(capfd):
    from multiprocessing import Queue

    capfd.readouterr()

    queue = Queue()

    for i in data:
        queue.put(i)
    queue.put(STOP)

    progress_bar_worker(queue, 5)

    std, err = capfd.readouterr()

    assert expected_bar.match(err)


def test_bar(capfd):
    capfd.readouterr()     # clean buffers

    with progress_bar(len(data)) as queue:
        for i in data:
            queue.put(i)

    assert queue.empty()

    # let's capture the output on both standard stream
    std, err = capfd.readouterr()

    # tqdm writes to err stream
    assert expected_bar.match(err)
