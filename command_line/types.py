from abc import ABC, abstractmethod
from argparse import ArgumentTypeError
from typing import Iterable, Any

from utils import abstract_property


class StringHandlingMixin(ABC):
    # TODO: this needs a documentation and some reasonable names

    @abstract_property
    def separator(self):
        pass

    @abstract_property
    def item_type(self):
        pass

    @abstract_property
    def data_type(self):
        pass

    def __init__(self, string):
        try:
            self.data = self.data_type(
                [
                    self.item_type(value) if value != '' else None
                    for value in string.split(self.separator)
                ]
                if self.separator else
                self.item_type(string)
            )
        except (TypeError, ValueError) as e:
            raise ArgumentTypeError(*e.args)


class Subset(ABC):

    @abstractmethod
    def get_iterator(self, iterable: Iterable[Any]) -> Iterable:
        return iterable

    def get(self, iterable: Iterable[Any]):
        return list(self.get_iterator(iterable))


def positive_int(value):
    value = int(value)
    if value < 0:
        raise ValueError('Indices need to be positive integers')
    return value


class Indices(Subset, StringHandlingMixin):

    separator = ','
    item_type = positive_int
    data_type = set

    def get_iterator(self, iterable):
        for i, value in enumerate(iterable):
            if i in self.data:
                yield value


class Slice(Subset, StringHandlingMixin):

    separator = ':'
    item_type = int
    data_type = tuple

    def get_iterator(self, iterable):
        return iterable[slice(*self.data)]


def one_of(*types):

    def one_of_types(string):

        exceptions = []
        for type_constructor in types:
            try:
                return type_constructor(string)
            except ArgumentTypeError as e:
                exceptions.append(f'{type_constructor.__name__}: {e}')

        names = ', '.join(t.__name__ for t in types)
        exceptions = ''.join('\n\t' + e for e in exceptions)

        raise ArgumentTypeError(
            f'Argument {string} does not match any of allowed types: {names}.\n' +
            f'Following exceptions has been raised: {exceptions}'
        )

    return one_of_types


def dsv(value_type, delimiter=','):
    """Delimiter Separated Values"""
    def closure(value):
        return [
            value_type(y)
            for y in value.split(delimiter)
        ]
    return closure
