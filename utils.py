from abc import ABCMeta, abstractmethod

import os

try:
    from numba import jit
except ImportError:
    def jit(func):
        return func

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable


class AbstractRegisteringType(ABCMeta):

    def __init__(cls, name, bases, attributes):
        super().__init__(name, bases, attributes)

        if not hasattr(cls, 'members'):
            cls.members = {}

        if hasattr(cls, 'name'):
            cls.members[cls.name] = cls
            for base in bases:
                if hasattr(base, 'name') and base.name in cls.members:
                    del cls.members[base.name]


def abstract_property(method):
    return property(abstractmethod(method))


def available_cores():
    # TODO test returns int
    return len(os.sched_getaffinity(0))
