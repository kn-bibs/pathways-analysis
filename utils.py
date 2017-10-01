from abc import ABCMeta, abstractmethod


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
