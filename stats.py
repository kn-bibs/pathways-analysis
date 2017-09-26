import numpy as np
from statsmodels.stats.weightstats import ttest_ind


def ttest_ind_mean(case: Phenotype, base_samples: Phenotype, alternative="two-sided"):
    """
    Two sided t-test of case sample(s) and mean expression values in base_samples across all genes
    :param case: data frame with case sample(s)
    :param base_samples: data frame with base sample(s)
    :param alternative: string with the alternative hypothesis, H1, has to be one of the following:
                        ‘two-sided’: H1: difference in means not equal to value (default)
                        ‘larger’ : H1: difference in means larger than value
                        ‘smaller’ : H1: difference in means smaller than value
    :return:
            tstat : float or numpy array in case of multiple case samples - test statisic
            pvalue : float or numpy array in case of multiple case samples - pvalue of the t-test
            df : int or float - degrees of freedom used in the t-test
    """
    l = [np.mean(row) for (idx, row) in base_samples.iterrows()]

    return ttest_ind(case, l, alternative)