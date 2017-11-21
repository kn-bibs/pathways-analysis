from collections import namedtuple
from contextlib import contextmanager
from multiprocessing import Manager, Queue, Process

from utils import tqdm
from utils import available_cores


def worker(func, input, progress_bar, output, *args):
    while True:
        data = input.get()
        # log że started

        if not data:
            continue

        if data == 'STOP':
            return

        func(data, *args)

        output.append(data)
        progress_bar.put(1)


def progress_bar_worker(queue, total):
    bar = tqdm(total=total)
    for step in iter(queue.get, None):
        if step == 'STOP':
            return
        if step:
            bar.update(step)


api_template = namedtuple('API', 'queue, results')


@contextmanager
def multiprocessing_queue(target, args, processes, total):
    m = Manager()
    results = m.list()
    queue = Queue()
    progress_queue = Queue()

    api = api_template(queue, results)

    progress = Process(target=progress_bar_worker, args=(progress_queue, total))

    processes_cnt = processes or available_cores()

    worker_args = [target, queue, progress_queue, results]
    if args:
        worker_args.extend(args)

    processes = [
        Process(target=worker, args=worker_args)
        for _ in range(processes_cnt)
    ]

    progress.start()

    yield api

    for process in processes:
        queue.put('STOP')

    for process in processes:
        process.start()

    for process in processes:
        process.join()

    progress_queue.put('STOP')
    progress.join()


class Pool:

    def __init__(self, processes):
        self.processes = processes

    def map(self, func, iterable, shared_args=tuple()):

        if self.processes == 1:
            # for profiling and debugging on process works better
            return map(lambda item: func(item, *shared_args), tqdm(iterable))

        with multiprocessing_queue(func, shared_args, self.processes, total=len(iterable)) as api:
            for item in iterable:
                api.queue.put(item)
        return api.results
