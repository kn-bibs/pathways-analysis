import operator

from multiprocess import Pool
from multiprocess import worker
from multiprocess import STOP


def test_imap():
    data = [1, 2, 3, 4]

    # for single process Pool has simplified logic,
    #  thus requires separate test iteration
    for num_processes in [1, 2]:

        pool = Pool(num_processes)

        # let's square: func = pow(x, 2)
        squared = pool.imap(pow, data, shared_args=[2])

        assert set(squared) == {1, 4, 9, 16}


def test_worker():
    from multiprocessing import Manager, Queue

    manager = Manager()
    output = manager.list()

    input_queue = Queue()
    progress_queue = Queue()

    for i in [0, 1, 2]:
        input_queue.put(i)
    input_queue.put(STOP)

    # let's add 5 to each number from the input queue
    worker(operator.add, input_queue, progress_queue, output, 5)

    assert list(output) == [5, 6, 7]
