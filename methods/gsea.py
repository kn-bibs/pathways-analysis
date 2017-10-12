from models import Experiment
from .method import Method


# TODO: do not use help, use docstring instead
class GSEA(Method):
    """Not finished yet."""

    name = 'gsea'

    help = """
    This method implementation is not finished yet :(
    
    Please refer & cite following publications:
        - Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550)
        - Mootha, Lindgren, et al. (2003, Nat Genet 34, 267-273)
    """

    def __init__(self):
        pass

    def run(self, experiment: Experiment):
        pass
