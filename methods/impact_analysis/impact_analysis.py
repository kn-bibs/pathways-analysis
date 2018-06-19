from databases import KEGGPathways
from methods.method import Method, MethodResult
from metrics import mean
from models import Experiment, Gene
from networkx import get_edge_attributes
from numpy import log2, isnan
from stats import ttest, hypergeom_distribution
from statsmodels.stats.multitest import multipletests
from .constants import *
import pandas as pd


class IAPathway:

    def __init__(self, data):
        self.IF = round(data['IF'], 3)
        self.pvalue = round(data['pvalue'], 3)
        self.name = data['name']
        self.FDR = round(data['FDR'], 3)
        self.Bonferroni = round(data['Bonferroni'], 3)


class ImpactAnalysisResult(MethodResult):
    columns = ['name', 'IF', 'pvalue', 'FDR', 'Bonferroni']


class ImpactAnalysis(Method):
    """
    Impact analysis is a pathways analysis method that considers the magnitude of each gene’s expression change,
    their type and position in the given pathways, and their interactions.

    Pipeline schema:
        1. Gene expression change between two conditions (logFC) is calculated.
        2. Differentially expressed genes are identified.
        3. KEGG Pathways database is searched for pathways containing at least one differentially expressed gene.
        4. An impact factor (IF) is calculated for each pathway incorporating parameters such as the normalized
           fold change of the differentially expressed genes, the statistical significance of the set of pathway genes,
           and the topology of the pathway.
        5. Multiple hypothesis testing adjustments of each pathway significance are performed by FDR and Bonferroni
           correction calculation.

    Additional method arguments allow specifying organism for database selection and the threshold for identification
    of differentially expressed genes.

    For more information, please refer to:
    Draghici S, Khatri P, Tarca AL, Amin K, Done A, et al. (2007),
    A systems biology approach for pathway level analysis. Genome Res 17: 1537–1545.

    """

    help = __doc__

    name = 'impact_analysis'

    def __init__(self, organism: str = 'Homo sapiens', threshold: float = 0.05, **kwargs):
        """

        Args:
            organism: organism name (ex. 'Homo sapiens', 'human')
            threshold: float: threshold for identification of differentially expressed genes
        """
        if threshold < 0 or threshold > 1:
            raise ValueError('Indices need to be in (0,1) range')
        self.threshold = threshold
        self.org = organism
        self.degs = None  # differentially expressed genes
        self.FC = None
        self.experiment_genes = None

    def run(self, experiment: Experiment) -> ImpactAnalysisResult:
        """

        Returns:
            list of pathways sorted by their impact factor. Each pathway in the list has values of FDR and
            Bonferroni corrections assigned.
        """

        # calculate fold change
        self.FC = experiment.calculate_fold_change()

        # remove genes for witch fold change cannot be calculated correctly
        experiment.exclude_genes(list(self.FC['FC'][isnan(self.FC['FC'])].index))

        # select differentialy expressed genes
        pvalue = ttest(experiment) <= self.threshold
        self.degs = pvalue[pvalue == True]

        if self.degs.size == 0:
            # if there are no DEGs anywhere, the problem of finding the impact on various pathways is meaningless
            print('No differentialy expressed genes.')
            return ImpactAnalysisResult([])

        self.experiment_genes = set([gene.name for gene in experiment.get_all().genes])

        db = KEGGPathways(self.org)
        pathways = {}

        for gene in [g.name for g in list(self.degs.index)]:
            ps = db.search_by_gene(gene)
            for (k, v) in ps.items():
                if k not in pathways.keys():
                    pathways[k] = v

        if not pathways:
            print('No pathways found in database.')
            return ImpactAnalysisResult([])

        res = pd.DataFrame(columns=['name', 'IF', 'pvalue'])
        for (code, descr) in pathways.items():
            pathway = db.get_pathway(code)
            impact_factor, pval = self.calculate_impact_factor(experiment, pathway)
            if impact_factor is not None and pval is not None:
                res.loc[len(res.index)] = [descr, impact_factor, pval]

        res['FDR'], res['Bonferroni'] = self.calculate_corrections(res['pvalue'])
        ifp_pathways = [IAPathway(res.loc[i]) for i in range(len(res.index))]
        ifp_pathways.sort(key=lambda x: x.IF if not isnan(x.IF) else 0, reverse=True)

        return ImpactAnalysisResult(ifp_pathways)

    def calculate_impact_factor(self, experiment, pathway):

        path_genes = set([x.strip() for x in ' ,'.join(pathway.nodes).split(',')])
        DEGs_set = set([gene.name for gene in list(self.degs.index)])

        # no DEGs in pathway
        if len(path_genes & DEGs_set) == 0:
            return None, None

        pval_path = hypergeom_distribution(len(path_genes & DEGs_set), len(self.experiment_genes), self.degs.size,
                                           len(path_genes & self.experiment_genes))

        if pval_path != 0:
            impact_factor = log2(pval_path)

            impact_factor += sum(
                [abs(self.calculate_perturbation_factor(experiment, gene, pathway)) for gene in pathway.nodes]) / len(
                path_genes & DEGs_set) * mean([abs(i) for i in self.FC['logFC'].values if not isnan(i)])

        else:
            impact_factor = MAX_IF

        return impact_factor, pval_path

    def calculate_perturbation_factor(self, experiment, gene, pathway, visited=None):

        visited = [] if not visited else visited

        if len(set(gene.split(',')) & set(self.experiment_genes)) == 0:
            return 0
        else:
            for name in gene.split(','):
                if name.strip() in self.experiment_genes:
                    # get ΔE (logFC)
                    pf = self.FC['logFC'][Gene(name)] if not isnan(self.FC['logFC'][Gene(name)]) else MAX_IF

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

    def calculate_corrections(self, pvalues):
        """

        Args:
            pvalues: array_like vector of p-values

        Returns:
            two arrays of p-values corrected for multiple tests using FDR and Bonferroni correction method

        """

        return multipletests(pvalues, method='fdr_bh')[1], multipletests(pvalues, method='bonferroni')[1]
