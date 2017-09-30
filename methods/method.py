from abc import abstractmethod
from models import Experiment
from utils import AbstractRegisteringType, abstract_property


class Method(metaclass=AbstractRegisteringType):

    def __init__(self, *args, **kwargs):
        pass

    @abstract_property
    def help(self):
        """Return string providing help for this method."""
        pass

    @abstract_property
    def name(self):
        """Return method name used internally and in command line interface."""
        pass

    @abstractmethod
    def run(self, experiment: Experiment):
        pass
