from abc import ABC, abstractmethod
from copy import copy
from random import shuffle

from methods.gsea.signatures import GeneSet
from models import SampleCollection, Experiment


def shuffle_and_divide(merged_collection, midpoint):
    shuffled = shuffle(merged_collection.samples)
    return SampleCollection(shuffled[:midpoint]), SampleCollection(shuffled[midpoint:])


class Shuffler(ABC):

    @abstractmethod
    def __init__(self, experiment: Experiment, rank, score):
        self.experiment = experiment
        self.rank = rank
        self.score = score
        self.gene_set = None

    def set_gene_set(self, gene_set: GeneSet):
        self.gene_set = gene_set

    @abstractmethod
    def permute_and_score(self):
        pass


class PhenotypeShuffler(Shuffler):

    def __init__(self, experiment: Experiment, rank, score):
        super().__init__(experiment, rank, score)
        self.all_samples = experiment.case + experiment.control
        self.cases_cnt = len(experiment.case.samples)

    def permute_and_score(self):
        random_case, random_control = shuffle_and_divide(self.all_samples, self.cases_cnt)
        ranked_list = self.rank(random_case, random_control)

        return self.score(ranked_list, self.gene_set)


class GeneShuffler(Shuffler):

    def __init__(self, experiment: Experiment, rank, score):
        super().__init__(experiment, rank, score)
        self.gene_labels = list(experiment.control.genes)
        self.permutation = copy(self.gene_labels)

    def permute_and_score(self):
        shuffle(self.permutation)

        ranked_list = self.rank(
            self.experiment.case, self.experiment.control,
            labels_map=dict(zip(self.gene_labels, self.permutation))
        )
        return self.score(ranked_list, self.gene_set)


