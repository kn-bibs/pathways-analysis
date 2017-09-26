from abc import ABCMeta, abstractmethod


class AbstractRegisteringType(ABCMeta):

    def __init__(cls, name, bases, attributes):
        super().__init__(name, bases, attributes)

        if not hasattr(cls, 'members'):
            cls.members = set()

        cls.members.add(cls)
        cls.members -= set(bases)


def abstract_property(method):
    return property(abstractmethod(method))
