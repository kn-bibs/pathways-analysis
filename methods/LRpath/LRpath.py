from methods.method import Method, MethodResult
from models import Experiment

class LRpathResult(MethodResult):
    columns = []
    description = """  """


class LRpath(Method):
    """
    method TODO
    """

    help = __doc__

    name = 'LRpath'

    legal_disclaimer = """  """


    def __init__(self):
        pass

    def run(self, experiment: Experiment) -> LRpathResult:
        return LRpathResult([])
