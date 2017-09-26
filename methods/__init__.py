from abc import ABC, abstractmethod
from models import Experiment


class Method(ABC):

    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def run(self, experiment: Experiment):
        pass
