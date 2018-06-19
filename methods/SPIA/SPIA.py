from methods.method import Method, MethodResult
from models import Experiment
from stats import ttest
from databases import KEGGPathways
from networkx import get_edge_attributes

#todo uzupełnić requirements

import numpy as np
import pandas as pd
from scipy import stats
import random
import math
from statsmodels.sandbox.stats import multicomp
from scipy.stats import norm

class SPIAPathway:

    def __init__(self, data):
        self.id = data['id']
        self.pNDE = round(data['pNDE'], 3)
        self.pPERT = round(data['pPERT'], 3)
        self.name = data['name']
        self.pG = round(data['pG'], 3)
        self.pGFWER = round(data['pGFWER'], 3)
        self.pGfdr = round(data['pGfdr'], 3)
        self.status = data['status']

class SPIAResult(MethodResult):
    #TODO description = """ """
    '''

    '''
    columns=['id','name','pNDE', 'pPERT', 'pG', 'pGfdr', 'pGFWER', 'status']

class SPIA(Method):
    #todo description!!!!
    '''

    '''
    help = __doc__
    name = "SPIA"

    def __init__(self, organism: str = 'hsa', threshold = 0.05):
        SPECIES = [["anopheles", "Anopheles gambiae", "Ag", "aga", "anoGam", "7165"],
                   ["bovine", "Bos taurus", "Bt", "bta", "bosTau", "9913"],
                   ["canine", "Canis familiaris", "Cf", "cfa", "canFam", "9615"],
                   ["chicken", "Gallus gallus", "Gg", "gga", "galGal", "9031"],
                   ["chimp", "Pan troglodytes", "Pt", "ptr", "PanTro", "9598"],
                   ["ecoliK12", "Escherichia coli K12", "EcK12", "eco", None, "562,83333,511145"],
                   ["ecoliSakai", "Escherichia coli Sakai", "EcSakai", "ecs", None, "83334"],
                   ["fly", "Drosophila melanogaster", "Dm", "dme", "dm", "7227"],
                   ["human", "Homo sapiens", "Hs", "hsa", "hg", "9606"],
                   ["mouse", "Mus musculus", "Mm", "mmu", "mm", "10090"],
                   ["pig", "Sus scrofa", "Ss", "ssc", "susScr", "9823"],
                   ["rat", "Rattus norvegicus", "Rn", "rno", "rn", "10116"],
                   ["rhesus", "Macaca mulatta", "Mmu", "mcc", "rheMac", "9544"],
                   ["worm", "Caenorhabditis elegans", "Ce", "cel", "ce", "6239"],
                   ["xenopus", "Xenopus laevis", "Xl", "xla", "NA", "8355"],
                   ["yeast", "Saccharomyces cerevisiae", "Sc", "sce", "sacCer", "4932,559292"],
                   ["zebrafish", "Danio rerio", "Dr", "dre", "danRer", "7955"]]

        for x in SPECIES:
            if organism in x:
                self.organism = x[1]
                break
        else:
            raise Exception("Unknown organism")

        #todo dodać inne zmienne do zmieniania
        if threshold < 0 or threshold > 1:
            raise ValueError('Indices need to be in (0,1) range')
        self.threshold = threshold


    @staticmethod
    def load_data_dict(dictionary,all):
        '''
        Part of this code is imported from https://github.com/iseekwonderful/PyPathway on MIT license.Load data from certain json.
        :param dictionary:  the json, this method will load data form json
        :param all: a python list of total genes. e.g. ['A', 'B', 'C', 'D']
        :return: the loaded data
        '''
        data=dictionary
        Bb=set(all)
        datpT = {}
        id2name = data['id2name']
        for pid, v in data.items():
            if pid == 'id2name':
                continue
            datpT[pid] = {}
            for key, d in v.items():
                if key == 'row_names':
                    for idname in range(len(d)):
                        Aaa=set(d[idname].split(','))
                        Ddd=list(Aaa.intersection(Bb))
                        if Ddd:
                            if len(Ddd)>1:
                                print("Error with gene names in pathway occured")
                            else:
                                d[idname]=Ddd[0]

                    datpT[pid][key] = d

                else:
                    m = np.zeros(
                        (len(v['row_names']), len(v['row_names'])))
                    for x in d:
                        m[x[0]][
                            x[1]] = 1
                    datpT[pid][
                        key] = m
        return datpT, id2name

    @staticmethod
    def calculate_spia(de, all, dictionary, nB=2000, beta=None, combine='fisher'):  # todo: opis
        '''
        This code is imported from https://github.com/iseekwonderful/PyPathway on MIT license. Contains SPIA alghoritm.
        Args:
            de: a python dict of DEGs. key:gene, value:fold-change. e.g. {'A':2.1, 'B':3.3 ...}
            all: a python list of total genes. e.g. ['A', 'B', 'C', 'D']
            dictionary: a json
            nB: number of iterations
            beta: list of gene relations values, if None the values are default
            combine:

        Returns:

        '''

        de = {k: float(v) for k, v in de.items()}
        all = [x for x in all]
        datpT_ALL, id2name = SPIA.load_data_dict(dictionary,all)

        # todo przerzucić to do osobnego pliku
        rel = ["activation", "compound", "binding/association", "expression", "inhibition",
               "activation_phosphorylation", "phosphorylation", "inhibition_phosphorylation",
               "inhibition_dephosphorylation", "dissociation", "dephosphorylation",
               "activation_dephosphorylation", "state change", "activation_indirect",
               "inhibition_ubiquination", "ubiquination", "expression_indirect",
               "inhibition_indirect", "repression", "dissociation_phosphorylation",
               "indirect_phosphorylation", "activation_binding/association",
               "indirect", "activation_compound", "activation_ubiquination"]
        inter_value = [1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1, -1, 0, 1, -1, -1, 0, 0, 1, 0, 1,
                       1] or beta
        rel_dict = {rel[i]: inter_value[i] for i in range(len(rel))}
        datp_ALL = {}
        for k, v in datpT_ALL.items():
            sizem = len(v[rel[0]][0])
            s, con = np.zeros((sizem, sizem)), np.zeros((sizem,sizem))
            for kk, vv in rel_dict.items():
                con += v[kk] * abs(vv)
                s += v[kk] * vv
            zz = np.reshape(np.repeat(con.sum(axis=0), sizem), (sizem, sizem))
            z = np.transpose(zz)
            z[z == 0] = -1
            r = np.divide(s, z)
            datp_ALL[k] = r
        smPFS, tAraw, tA, pNDE, pb, pG, status = {}, {}, {}, {}, {}, {}, {}
        # calculate the Ac
        for k, v in datp_ALL.items():
            row_names = datpT_ALL[k]['row_names']
            # let first calculate the pNDE
            noMy = len(
                set(row_names) & set(de.keys()))
            pNDE[k] = stats.hypergeom.sf(noMy - 1, len(all), len(set(row_names) & set(all)), len(de))
            # then calculate the Ac and pPERT
            M = np.eye(v.shape[0]) * -1 + v
            if np.linalg.det(M) == 0:
                smPFS[k], tAraw[k], tA[k], pb[k] = np.nan, np.nan, np.nan, np.nan
                continue
            X = []
            for x in row_names:
                if x in de:
                    X.append(de[x])
                else:
                    X.append(0)
            pfs = np.linalg.solve(M, -np.array(X))
            smPFS[k] = sum(pfs - X)
            tAraw[k] = smPFS[k]
            pfstmp = []
            de_sample = list(de.values())
            all_sample = [i for i, x in enumerate(row_names) if x in all]
            length = len(X)
            for i in range(nB):
                x = np.zeros(length)
                sp = random.sample(de_sample, noMy)
                idx = random.sample(all_sample, noMy)
                x[idx] = sp
                tt = np.linalg.solve(M, -x)
                pfstmp.append(sum(tt - x))
            tA[k] = tAraw[k] - np.median(np.array(pfstmp))
            if tA[k] > 0:
                status[k] = "Activated"
            else:
                status[k] = "Inhibited"
            ob = tA[k]
            pfstmp = np.array(pfstmp) - np.median(np.array(pfstmp))
            if ob > 0:
                pb[k] = sum([1 for pf in pfstmp if pf >= ob]) / len(pfstmp) * 2
                if pb[k] <= 0:
                    pb[k] = 1 / nB / 100
                elif pb[k] >= 1:
                    pb[k] = 1
            elif ob < 0:
                pb[k] = sum([1 for pf in pfstmp if pf <= ob]) / len(pfstmp) * 2
                if pb[k] <= 0:
                    pb[k] = 1 / nB / 100
                elif pb[k] >= 1:
                    pb[k] = 1
            else:
                pb[k] = 1
            if combine == 'fisher':
                c = pNDE[k] * pb[k]
                pG[k] = c - c * math.log(c)
            else:
                pG[k] = norm.cdf(norm.ppf(pNDE[k]) + norm.ppf(pb[k]) / math.sqrt(2))
        _, o, _, _ = multicomp.multipletests(list(pG.values()), method='fdr_bh')
        pGfdr = {list(pG.keys())[i]: o[i] for i in range(len(list(pG.keys())))}
        _, o, _, _ = multicomp.multipletests(list(pG.values()), method='bonferroni')
        pGbf = {list(pG.keys())[i]: o[i] for i in range(len(list(pG.keys())))}
        df=pd.DataFrame(columns=['id','name','pNDE', 'pPERT', 'pG', 'pGfdr', 'pGFWER', 'status'])
        for k in pG:
            df.loc[len(df.index)] = [k,id2name[k],pNDE[k],pb[k],pG[k],pGfdr[k],pGbf[k],status[k]]
        df2=[SPIAPathway(df.loc[i]) for i in range(len(df.index))]
        return df2

    def run(self, experiment: Experiment) -> SPIAResult:

        pvalue = ttest(experiment) <= self.threshold
        calc_f = experiment.calculate_fold_change()
        all=pvalue.index.tolist()
        for a in range(len(all)):
            all[a]=all[a].name
        p= pvalue[pvalue == True].index.tolist()
        de={}
        for i in p:
            de[i.name] = calc_f["FC"][i]
        json={}
        if len(de) == 0:
            # if there are no DEGs anywhere, the problem of finding the impact on various pathways is meaningless
            print('No differentialy expressed genes.')
            return SPIAResult([])
        db = KEGGPathways(self.organism)
        pathways = {}
        for gene in de.keys():
            ps = db.search_by_gene(gene)
            for (k, v) in ps.items():
                if k not in pathways.keys():
                    pathways[k] = v
        if not pathways:
            print('No pathways found in database.')
            return SPIAResult([])
        for (id, descr) in pathways.items():
            pathway = db.get_pathway(id)
            path_genes=set(pathway.nodes)
            path_genes=list(path_genes)

            #todo: to też zmienić
            rel = ["activation", "compound", "binding/association", "expression", "inhibition",
                   "activation_phosphorylation", "phosphorylation", "inhibition_phosphorylation",
                   "inhibition_dephosphorylation", "dissociation", "dephosphorylation",
                   "activation_dephosphorylation", "state change", "activation_indirect",
                   "inhibition_ubiquination", "ubiquination", "expression_indirect",
                   "inhibition_indirect", "repression", "dissociation_phosphorylation",
                   "indirect_phosphorylation", "activation_binding/association",
                   "indirect", "activation_compound", "activation_ubiquination"]

            #todo beta jeśli nieznane interakcje

            interaction_list={i:[] for i in rel}
            x=get_edge_attributes(pathway, 'type')
            for gene1,interaction in x.items():
                interaction = '_'.join(interaction)
                if interaction in interaction_list.keys():
                    interaction_list[interaction].append([path_genes.index(gene1[0]),path_genes.index(gene1[1])])
                else:
                    interaction_list[interaction]=[[path_genes.index(gene1[0]), path_genes.index(gene1[1])]]
            interaction_list['row_names']=path_genes
            json[id]=interaction_list
        json['id2name'] = pathways
        s=SPIA.calculate_spia(de,all,json)
        return SPIAResult(s)