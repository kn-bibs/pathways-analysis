import random
from itertools import chain
from math import sqrt
from operator import itemgetter
from textwrap import dedent
from typing import List, Iterable

import numpy as np

from declarative_parser.parser import action, Argument
from declarative_parser.types import positive_int
from utils import jit

import multiprocess
from methods.gsea.shufflers import PhenotypeShuffler, GeneShuffler
from methods.method import Method, MethodResult
from models import Experiment, SampleCollection
from .metrics import SignalToNoise, ranking_metrics
from .signatures import DatabaseParser, GeneSet


class GSEAResult(MethodResult):

    # TODO: FWER p-values
    columns = ['name', 'enrichment', 'nominal_p_value', 'fdr']

    description = """
    For help with interpreting the data see:
    http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_GSEA_Results
    """


class ScoreDistribution:

    __slots__ = ('negative_scores', 'positive_scores')

    def __init__(self, scores: Iterable=None):
        self.negative_scores = []
        self.positive_scores = []
        if scores:
            for score in scores:
                self.append(score)

    def append(self, score):
        if score > 0:
            self.positive_scores.append(score)
        elif score < 0:
            self.negative_scores.append(score)
        else:
            # TODO: is it the best way to do this?
            target = random.choice((self.positive_scores, self.negative_scores))
            target.append(score)

    def __iter__(self):
        return chain(self.negative_scores, self.positive_scores).__iter__()

    def __len__(self):
        return len(self.positive_scores) + len(self.negative_scores)


@jit
def is_more_extreme(x, enrichment):
    """Is x more extreme (more negative or more positive) than provided enrichment?"""
    return abs(x) > abs(enrichment)


class JavaGSEA(Method):
    # TODO: create wrapper for Desktop version from Broad Institute
    pass


class GeneralisedGSEA(Method):
    """
    GSEA method evaluates list of provided gene sets,
    looking for such gene sets which are significantly
    enriched in first of provided collections of samples (or class),
    when compared to the second "control" group of samples.

    The result is presented as list of gene sets sorted by their
    NES (normalized enrichment score) values;
    Additionally FDR q-values and nominal p-values
    are presented for each set of genes.
    (TODO: FWER: p-values)

    Schematic of pipeline:
        1. Ranked list of genes present in provided samples is created,
           with the ranking based on differential expression (using
           user-chosen ranking metric, default SignalToNoise).
        2. For each set of genes:
            - enrichment score is calculated,
              using KS-test-like walking sum statistic;
            - null (random) distribution of enrichment scores is created
              using gene or class permutations
            - nominal p-value is assessed using the null distribution
            - normalized enrichment score is calculated
        3. Multiple hypothesis testing adjustment is performed by FDR
           calculation, exploring ratios of values more extreme than
           normalized enrichment scores in the observed and
           in all random, normalized distributions.

    Required arguments:
        One has to choose database with gene sets to be used for analysis

    Additional arguments include:
     - permutations type and count adjustment
     - control of "p" exponent used in enrichment weighting for ranked lists creation
     - multiprocessing control
     - ranking metric

    Please refer & cite following publications:
        - Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550)
        - Mootha, Lindgren, et al. (2003, Nat Genet 34, 267-273)

    Please use --show_licence to display licence and copyright details.
    """

    help = __doc__

    name = 'gsea'

    legal_disclaimer = """
        Molecular signature databases (MSigDB) are protected by 
        copyright (c) 2004-2017 Broad Institute, Inc., Massachusetts
        Institute of Technology, and Regents of the University of
        California, subject to the terms and conditions of the
        Creative Commons Attribution 4.0 International License:
        https://creativecommons.org/licenses/by/4.0/
    """

    database = DatabaseParser()

    shufflers = {'phenotypes': PhenotypeShuffler, 'genes': GeneShuffler}

    ranking_metric = Argument(
        type=lambda name: ranking_metrics[name],
        choices=list(ranking_metrics.keys()),
        default=SignalToNoise
    )

    # TODO: metrics
    def __init__(self, database, ranked_list_weight: float=1, ranking_metric=SignalToNoise,
                 permutation_type='genes', normalize_es=True,
                 processes: positive_int=0, permutations: positive_int=1000, min_genes=15,
                 max_genes=500, descending_sort=True, **kwargs):
        """

        Args:
            database: Database or object with database property.
            ranked_list_weight:
                an enrichment weighting exponent (p in publication)
                to control weight
                of producing ranked list, default 1
            permutation_type: 'genes' or 'phenotypes' (one of shufflers)
            permutations: count of permutation rounds to execute (the more, the mor accurate is pvalue estimation)
            normalize_es: should enrichment scores be normalized (adjusting for variation in gene sets)?
            processes: a number of processes to use; by default all available cores will be utilized

        """
        # TODO: store opts in a dict?
        if hasattr(database, 'database'):
            database = database.database
        self.database = database
        self.ranked_list_weight = ranked_list_weight
        self.calculate_rank = ranking_metric()
        self.normalize_es = normalize_es
        self.processes = processes
        self.permutations = permutations
        self.gene_sets = [s for s in self.database.gene_sets.values() if min_genes > len(s.genes) < max_genes]
        self.descending_sort = descending_sort

        self.shuffler_class = self.shufflers[permutation_type]
        self.shuffler = None
        # TODO: p-value, fdr cutoff

    def run(self, experiment: Experiment):
        """Return list of gene sets sorted by normalized enrichment score.

        Each gene set in the list has FDR and enrichment_score assigned."""

        # L in the Subramanian2005 publication
        ranked_list = self.create_ranked_gene_list(
            experiment.case, experiment.control
        )
        # gene_set is S in the publication
        self.shuffler = self.shuffler_class(
            experiment,
            self.create_ranked_gene_list,
            self.calculate_enrichment_score,
        )

        args = (ranked_list, )

        pool = multiprocess.Pool(self.processes)
        gene_sets = pool.imap(self.analyze_gene_set, self.gene_sets, shared_args=args)

        sorted_gene_sets = sorted(gene_sets)

        self.compute_fdr(sorted_gene_sets)

        return GSEAResult(sorted_gene_sets)

    def analyze_gene_set(self, gene_set: GeneSet, ranked_list):
        # 1. step in the publication (Calculation of an Enrichment Score)
        enrichment_score = self.calculate_enrichment_score(ranked_list, gene_set)

        # 2. step in the publication (Estimation of Significance Level of ES)
        null_distribution = self.enrichments_for_permuted_labels(gene_set)
        # significance level is a nominal p-value here
        nominal_p_value = self.estimate_significance_level(enrichment_score, null_distribution)

        # 3. step in the publication (Adjustment for Multiple Hypothesis Testing)
        normalized_enrichment, normalized_null_distribution = self.normalize_enrichment(
            enrichment_score,
            null_distribution
        )

        # warning: FDR will be calculated once all gene sets are analysed!

        gene_set.enrichment = normalized_enrichment
        gene_set.nominal_p_value = nominal_p_value
        gene_set.null_distribution = null_distribution

        return gene_set

    @action
    def show_licence(namespace):
        """Print out licence and legal information."""
        print(dedent(GeneralisedGSEA.legal_disclaimer))

    @action
    def download_all(namespace):
        """Fetch all databases"""
        pass

    def create_ranked_gene_list(self, case: SampleCollection, control: SampleCollection, labels_map=None):
        """

        Args:
            case:
            control:
            labels_map: used for permutations or if sample collection does not have labels

        Returns:

        """
        assert case.genes == control.genes

        genes = case.genes

        if not labels_map:
            def map_label(gene):
                return gene
        else:
            def map_label(gene):
                return labels_map[gene]

        # case/ control should be shuffled, not genes! this won't work for continuous

        # dict in 3.6 are ordered!
        return sorted(
            [
                (map_label(gene), self.calculate_rank(case.of_gene(gene), control.of_gene(gene)))
                for gene in genes
            ],
            key=itemgetter(1),
            reverse=self.descending_sort
        )

    @jit
    def calculate_enrichment_score(self, ranked_list, gene_set: GeneSet):
        # TODO: review of formulas more than welcome;
        # based on formulas from "Appendix: Mathematical Description of Methods"
        maximum_deviation = 0
        running_sum_statistic_hits = 0
        running_sum_statistic_misses = 0

        p = self.ranked_list_weight

        n = len(ranked_list)
        nh = len(gene_set)

        # P_miss(S, i, j)
        decrement = 1 / (n - nh)

        # weight, N_R
        hit_denominator = 0
        for gene, rank in ranked_list:
            if gene in gene_set:
                hit_denominator += abs(pow(rank, p))

        for gene, rank in ranked_list:

            power_of_rank = pow(rank, p)

            # hit
            if gene in gene_set:
                increment = power_of_rank / hit_denominator
                running_sum_statistic_hits += increment
            # miss
            else:
                running_sum_statistic_misses -= decrement

            diff = running_sum_statistic_hits - running_sum_statistic_misses

            if abs(diff) > abs(maximum_deviation):
                maximum_deviation = diff

        return maximum_deviation

    def enrichments_for_permuted_labels(self, gene_set):
        """Create null distribution by repetitive permutations of gene labels"""
        # es_null is used to store a histogram representing null distribution
        es_null = ScoreDistribution()

        shuffler = self.shuffler
        shuffler.set_gene_set(gene_set)

        for score in range(self.permutations):
            score = shuffler.permute_and_score()
            es_null.append(score)

        return es_null

    @staticmethod
    @jit
    def estimate_significance_level(enrichment_score, null_distribution: ScoreDistribution):
        """Estimate nominal p-value using only this side (tail) of null (random) distribution

        that corresponds to the sign of provided enrichment score.

        Value 0 indicates that the nominal p-value is < (1 / self.permutations).
        """
        if enrichment_score > 0:
            null = null_distribution.positive_scores
        else:
            null = null_distribution.negative_scores

        abs_score = abs(enrichment_score)

        hits = 0

        for random_score in null:
            if abs(random_score) > abs_score:
                hits += 1

        if not null:
            return 0

        p_value = hits / len(null)

        return p_value

    @staticmethod
    def compute_fdr(analyzed_gene_sets: List[GeneSet]):
        """FDR = a ratio of more extreme results in random distributions / observed.

        Citing from GSEA documentation:

            The FDR is a ratio of two distributions:
            (1) the actual enrichment score versus the enrichment scores
                for all gene sets against all permutations of the dataset and
            (2) the actual enrichment score versus the enrichment scores
                of all gene sets against the actual dataset.

        For formal definition, see last paragraph of:
        "Appendix: Mathematical Description of Methods".

        Args:
            analyzed_gene_sets: already analyzed gene sets
        """

        for gene_set in analyzed_gene_sets:

            normalized_enrichment = gene_set.enrichment

            more_extreme_random = 0
            all_random = 0

            more_extreme_observed = 0
            observations = 0

            for other_set in analyzed_gene_sets:
                if other_set == gene_set:
                    continue

                more_extreme_random += sum(
                    1 for random_enrichment in other_set.null_distribution
                    if is_more_extreme(random_enrichment, normalized_enrichment)
                )
                all_random += len(other_set.null_distribution)

                if is_more_extreme(other_set.enrichment, normalized_enrichment):
                    more_extreme_observed += 1
                observations += 1

            # this controls division by zero and provides a shortcut to quit if there are no results
            if not more_extreme_random:
                # TODO: less than
                gene_set.fdr = 0
                continue

            nominator = more_extreme_random / all_random
            denominator = more_extreme_observed / observations

            gene_set.fdr = nominator / denominator if denominator else None

    def normalize_enrichment(self, enrichment_score: float, null_distribution: ScoreDistribution):
        """Normalize enrichment by dividing be mean.

        See: Multiple Hypothesis Testing, point 3 in publication.
        """
        if not self.normalize_es:
            return null_distribution, enrichment_score

        null_positive = null_distribution.positive_scores
        null_negative = null_distribution.negative_scores

        mean_positive = np.mean(null_positive) if null_positive else None
        mean_negative = np.mean(null_negative) if null_negative else None

        mean_negative_abs = abs(mean_negative) if mean_negative else None

        def normalized(score):
            if score > 0:
                score /= mean_positive
            else:
                score /= mean_negative_abs
            return score

        normalized_null = ScoreDistribution(
            (normalized(raw_score) for raw_score in null_distribution)
        )

        # TODO: what if mean_negative_abs is None but score is negative?
        # or mean positive is None and score positive?
        # ODP: GSEA 3.0 just stipulates NaN or blank
        # http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/FAQ#What_does_it_mean_for_a_gene_set_to_have_NES_and_nominal_p-values_of_NaN_.28also_shown_as_blanks.29.3F

        return normalized(enrichment_score), normalized_null


class SimpleGSEA(GeneralisedGSEA):
    """Implementation of GSEA as described in:

    Mootha, V. K., Lindgren, C. M., Eriksson, K. F., Subramanian, A., Sihag, S., Lehar, J.,
    Puigserver, P., Carlsson, E., Ridderstrale, M., Laurila, E., et al. (2003) Nat. Genet. 34,
    267–273.

    NOT TESTED
    """

    help = __doc__
    name = 'simple_gsea'

    database = DatabaseParser()

    # TODO: test this
    def calculate_enrichment_score(self, ranked_list, gene_set):
        # variable names were chosen to reflect description Supporting Text of GeneralisedGSEA

        maximum_deviation = 0
        running_sum_statistic = 0

        n = len(ranked_list)
        nh = len(gene_set)

        increment = sqrt(n - nh) / nh
        decrement = sqrt(nh / (n - nh))

        for gene, _ in ranked_list:
            if gene in gene_set:
                running_sum_statistic += increment
            else:
                running_sum_statistic -= decrement
            if abs(running_sum_statistic) > abs(maximum_deviation):
                maximum_deviation = running_sum_statistic

        return maximum_deviation

