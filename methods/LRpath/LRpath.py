import math
import numpy
import os.path
import sys

import pandas as pd
import statsmodels.api as sm
from scipy import stats

from databases import KEGGPathways
from methods.method import Method, MethodResult
from models import Experiment


class LRpathResult(MethodResult):
    columns = ['Path_ID', 'nGene', 'LRcoeff', 'odds_ratio', 'LRpvalue', 'catsigIDs']
    description = """ 
    For more information see:
    http://lrpath.ncibi.org/method.pdf 
    """


class LRPathway:

    def __init__(
            self, Path_ID: str, nGene: int, LRcoeff: float, odds_ratio: float, LRpvalue: float, catsigIDs: str
    ):
        self.LRcoeff = LRcoeff
        self.odds_ratio = odds_ratio
        self.LRpvalue = LRpvalue
        self.nGene = nGene
        self.Path_ID = Path_ID
        self.catsigIDs = catsigIDs


class LRpath(Method):
    """
    LRpath performs gene set enrichment testing, an approach used to test for predefined
    biologically-relevant gene sets that contain more significant genes from an experimental
    dataset than expected by chance. Given a high-throughput dataset with continuous significance
    values (i.e. p-values), LRpath tests for gene sets (termed concepts) that have significantly
    higher significance values (e.g. for differential expression) than expected at random.  Genes
    are mapped to concepts using their Entrez Gene IDs. The use of logistic regression allows the
    data to remain on a continuous scale while maintaining the interpretation of results in terms
    of an odds ratio , as is used with the standard Fisher's Exact test.

    Additional arguments include:
     - minimum number of unique gene IDs analyzed in category to be tested
     - maximum number of unique gene IDs analyzed in category to be tested
     - cutoff = Entrez gene IDs in each category with p-values < will be tested
     - database to be tested
     - lower and upper p-values to be used

    Please refer & cite following publications:
        - 'LRpath: a logistic regression approach for identifying enriched biological groups in gene expression data.'
            Sartor, Leikauf, Medvedovic (2009, Bioinformatics 25, 211-217)

    Please use --show_licence to display licence and copyright details.
    """

    help = __doc__

    name = 'LRpath'

    legal_disclaimer = """ Copyright 2010 The University of Michigan  """

    def __init__(
            self, database, min_g=10, max_g=None,
            cutoff=0.05, odds_min=0.001, odds_max=0.5,
    ):
        """

        Args:
            database: file with columns separated by tab, [1] first column should has pathway_id in first place (it can
                some others infromacions separated by spaces), and [2] second column should has gene exist in this
                pathway separated by space.
            min_g: the minimum number of unique gene IDs analyzed in category to be tested
            max_g: the maximum number of unique gene IDs analyzed in category to be tested
            cutoff: entrez gene IDs in each category with p-values < cutoff, will be tested
            odds_min: lower p-values be used
            odds_max: upper p-values to be used
        """

        self.database = database
        self.min_g = min_g
        self.max_g = max_g
        self.cutoff = cutoff
        self.odds_min = odds_min
        self.odds_max = odds_max

    def run(self, experiment: Experiment) -> LRpathResult:
        """
        Return list of results

        """

        geneids = experiment.case.genes
        match = experiment.case.as_array()
        data, names_sample = self.create_data(match)
        db = self.create_database()
        data, geneid = self.name_geneid(data, geneids)
        results = self.calc_siggenes(data, names_sample, geneid, db)

        return LRpathResult(results)

    @staticmethod
    def create_data(match):
        """"
        Create sample data and list of sample names

        """
        diki = {}
        names = []

        for name in match:
            names.append(name)
            col = []
            for i in match.index:
                test = (float(match.loc[[i]][name]))
                col.append(float(stats.chi2.pdf(test, 1)))
            diki[name] = col

        data = pd.DataFrame(diki)
        data.index = match.index

        return data, names

    def create_database(self):
        """
        Open file with database and initialize creating database accepted by method

        """
        if os.path.exists(self.database):
            if type(self.database) is not dict:
                db = self.get_list_db()

            else:
                db = self.database
            return db

        else:
            sys.exit("Path with database file is not exist")

    @staticmethod
    def name_geneid(data, geneids):
        geneid = []

        for geny in geneids:
            geneid.append(KEGGPathways().get_gene_code(gen=geny.name).split()[0])

        data['gene_name'] = data.index
        data.index = geneid

        return data, geneid

    def get_list_db(self):
        """
        Returns:
            A database accepted in LRpath method

        """
        file = open(self.database).readlines()
        db = {}

        for line in file:
            line = line.split('\t')
            i_name = line[0].split()[0]
            i_list = line[1].split()
            db[i_name] = i_list

        return db

    def calc_siggenes(self, data, names, geneid, database):
        """
        Calculate logistic regression from data.
         Args:
             data: DataFrame
             names: list of simple names
             geneid: list of geneid
             database: dick

        """
        name = names[0]
        signes = []
        for si in data[name]:
            if si == 0:
                si = 10 ** (-15)
            elif si == float('inf'):
                si = 1
            signes.append(si)
        data.update(pd.Series(signes, name=name, index=data.index))

        uniqids = list(set(geneid))
        numuniq = len(uniqids)
        lor_mult = math.log(self.odds_min) - math.log(self.odds_max)
        data['nlp'] = None

        for sig in data.index:
            up = pd.Series([float((-1) * math.log(float(data.loc[[sig]][name])))], name='nlp', index=[sig])
            data.update(up)

        newp = [None] * numuniq

        for num in range(numuniq):
            current = []
            for genid in geneid:
                if genid == uniqids[num]:
                    current.append(data.loc[[genid]]['nlp'][0])
            for cur in current:
                if cur is not None:
                    newp[num] = numpy.mean(cur)

        catsizes = {}
        yy = {}
        for k, v in database.items():
            catsizes[k] = len(v)
            if catsizes[k] >= self.min_g:
                yy[k] = v

        siggenes = []
        for i in range(len(uniqids)):
            if newp[i] is not None:
                if math.exp(-newp[i]) < self.cutoff:
                    siggenes.append(uniqids[i][4:])

        if self.max_g is None:
            self.max_g = 99999

        ind = 0
        catsigIDs = {}
        results = []
        for key, value in yy.items():
            path = []
            for el in value:
                for u in range(len(uniqids)):
                    if el == uniqids[u][4:]:
                        path.append(u)

            length = len(path)
            if self.max_g >= length >= self.min_g:
                ind = ind + 1
                cat = [0] * numuniq
                for s in path:
                    cat[s] = 1

                logmod = sm.GLM(pd.DataFrame(cat), pd.DataFrame(newp), family=sm.families.Binomial())
                t = logmod.fit()

                id_list = set()

                for gen in value:
                    if gen in siggenes:
                        id_list.add(int(gen))
                catsigIDs[ind] = ', '.join(map(str, id_list))

                if catsigIDs[ind]:
                    lrpath = LRPathway(
                        key, length,
                        float(t.params),
                        float(math.exp(lor_mult * float(t.params))),
                        float(t.pvalues),
                        catsigIDs[ind])
                    results.append(lrpath)

        return results
