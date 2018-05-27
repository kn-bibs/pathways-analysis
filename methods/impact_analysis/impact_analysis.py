from methods.method import Method, MethodResult
from models import Experiment


class ImpactAnalysisResult(MethodResult):
    columns = []
    description = """ test """


class ImpactAnalysis(Method):
    """
    method and workflow description and sources
    """

    help = __doc__

    name = 'impact_analysis'

    legal_disclaimer = """ test """

    def __init__(self):
        pass

    def run(self, experiment: Experiment) -> ImpactAnalysisResult:
        return ImpactAnalysisResult([])
