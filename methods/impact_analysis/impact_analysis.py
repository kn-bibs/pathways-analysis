from databases import KEGGPathways
from methods.method import Method, MethodResult
from metrics import mean
from models import Experiment, Gene
from networkx import get_edge_attributes
from numpy import log2, isnan, nan
from stats import ttest, hypergeom_distribution
from .constants import *


class IAPathway:

    def __init__(self, IF: float, pvalue: float, name: str):
        self.IF = IF
        self.pvalue = pvalue
        self.name = name

    def set_corrections(self, fdr: float, bonf: float):
        self.FDR = fdr
        self.Bonferroni = bonf


class ImpactAnalysisResult(MethodResult):
    # columns = ['name', 'IF', 'pvalue', 'FDR', 'Bonferroni'] #TODO
    columns = ['name', 'IF', 'pvalue']
    description = """ test """  # TODO


class ImpactAnalysis(Method):
    """
    method and workflow description and sources
    """  # TODO

    help = __doc__

    name = 'impact_analysis'

    legal_disclaimer = """ test """  # TODO

    def __init__(self, organism: str = 'Homo sapiens', threshold: float = 0.05, **kwargs):
        if threshold < 0 or threshold > 1:
            raise ValueError('Indices need to be in (0,1) range')
        self.threshold = threshold
        self.org = organism
        self.DEGs = None
        self.FC = None
        self.experiment_genes = None

    def run(self, experiment: Experiment) -> ImpactAnalysisResult:

        # select differentialy expressed genes
        pvalue = ttest(experiment) <= self.threshold
        self.DEGs = pvalue[pvalue == True]

        if self.DEGs.size == 0:
            # if there are no DEGs anywhere, the problem of finding the impact on various pathways is meaningless
            print('No differentialy expressed genes.')
            return ImpactAnalysisResult([])

        # calculate fold change
        self.FC = experiment.calculate_fold_change()

        self.experiment_genes = set([gene.name for gene in experiment.get_all().genes])

        db = KEGGPathways(self.org)
        pathways = {}

        for gene in [g.name for g in list(self.DEGs.index)]:
            ps = db.search_by_gene(gene)
            for (k, v) in ps.items():
                if k not in pathways.keys():
                    pathways[k] = v

        if not pathways:
            print('No pathways found in database.')
            return ImpactAnalysisResult([])

        ifp_pathways = []
        for (code, descr) in pathways.items():
            pathway = db.get_pathway(code)
            impact_factor, pval = self.calculate_impact_factor(experiment, pathway)
            if impact_factor is not None and pval is not None:
                ifp = IAPathway(impact_factor, pval, descr)
                ifp_pathways.append(ifp)

        ifp_pathways.sort(key=lambda x: x.IF if not isnan(x.IF) else 0, reverse=True)
        return ImpactAnalysisResult(ifp_pathways)

    def calculate_impact_factor(self, experiment, pathway):

        path_genes = set([x.strip() for x in ' ,'.join(pathway.nodes).split(',')])
        DEG_genes = set([gene.name for gene in list(self.DEGs.index)])

        # no DEGs in pathway
        if len(path_genes & DEG_genes) == 0:
            return None, None

        pval_path = hypergeom_distribution(len(path_genes & DEG_genes), len(self.experiment_genes), self.DEGs.size,
                                           len(path_genes & self.experiment_genes))

        if pval_path != 0:
            impact_factor = log2(pval_path)

            impact_factor += sum(
                [self.calculate_perturbation_factor(experiment, gene, pathway) for gene in pathway.nodes]) / len(
                path_genes & DEG_genes) * mean([abs(i) for i in self.FC['FC'].values if not isnan(i)])

        else:
            impact_factor = MAX_IF

        return round(impact_factor, 3), round(pval_path, 3)

    def calculate_perturbation_factor(self, experiment, gene, pathway, visited=None):

        visited = [] if not visited else visited

        if len(set(gene.split(',')) & set(self.experiment_genes)) == 0:
            return 0
        else:
            for name in gene.split(','):
                if name.strip() in self.experiment_genes:
                    # get ΔE if measurable #TODO diffrent measures of expression change
                    pf = self.FC['FC'][Gene(name)] if not isnan(self.FC['FC'][Gene(name)]) else 0
                    # genes directly upstream
                    for edge in pathway.in_edges(gene):
                        if edge[0] not in visited:
                            beta = mean([interaction_weights[t] if t in interaction_weights.keys() else 0 for t in
                                         get_edge_attributes(pathway, 'type')[edge]])
                            # genes directly downstream
                            dstream = len(pathway.out_edges(edge[0]))
                            pf += self.calculate_perturbation_factor(experiment, edge[0], pathway,
                                                                     visited + [edge[0]]) * beta / dstream
                    return pf
