from models import Experiment
import pandas as pd
from scipy.stats import ttest_ind, hypergeom


def ttest(experiment: Experiment):
    """This is a two-sided test for the null hypothesis that 2 independent samples have identical average (expected) values.
        This test assumes that the populations have identical variances by default.

    Args:
        experiment: Experiment object with case and control sample(s)

    Returns:
        one-dimensional labeled array of p-values with Gene objects as labels

    """
    pvals = {}
    for (idx, row) in experiment.get_all().as_array().iterrows():
        control = [row[label] for label in experiment.control.labels]
        case = [row[label] for label in experiment.case.labels]
        pvals[idx] = ttest_ind(control, case).pvalue
    return pd.Series(pvals, name='p-value')


def hypergeom_distribution(k: int, m: int, n: int, s: int):
    """
    Survival function of hypergeometric distribution.

    Returns:
        Probability of drawing at least 'k' objects of Type I from bin,
        containing 'm' objects (of which 'n' are Type I objects),
        while 's' objects are randomly drawn without replacement from the total population.

    """
    h = hypergeom.sf(k - 1, m, n, s)
    return h
