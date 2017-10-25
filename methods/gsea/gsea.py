from textwrap import dedent

from command_line.parser import action
from models import Experiment
from methods.method import Method

from .signatures import DatabaseParser


# TODO: do not use help, use docstring instead
class GSEA(Method):
    """Not finished yet."""

    name = 'gsea'

    help = """
    Expected results, based on abstract
    
    Schematic of pipeline
    
    Required args description
    
    List of additional args?
    
    Please refer & cite following publications:
        - Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550)
        - Mootha, Lindgren, et al. (2003, Nat Genet 34, 267-273)
        
    Please use --show_licence to display licence and copyright details.
    """

    legal_disclaimer = """
        Molecular signature databases (MSigDB) are protected by 
        copyright (c) 2004-2017 Broad Institute, Inc., Massachusetts
        Institute of Technology, and Regents of the University of
        California, subject to the terms and conditions of the
        Creative Commons Attribution 4.0 International License:
        https://creativecommons.org/licenses/by/4.0/
    """

    database = DatabaseParser()

    def __init__(self, database, gene_set, **kwargs):
        """

        Args:
            database: Database or object with database property.
        """
        if hasattr(database, 'database'):
            database = database.database
        self.gene_set = database.gene_sets[gene_set]

    def run(self, experiment: Experiment):
        print(self.gene_set)

    @action
    def show_licence(namespace):
        """Print out licence and legal information."""
        print(dedent(GSEA.legal_disclaimer))

    @action
    def download_all(namespace):
        """Fetch all databases"""
        pass

    def calculate_rank(self):
        #
        pass

