from models import SampleCollection, Sample
import numpy as np
from statsmodels.stats.weightstats import ttest_ind
from typing import Union


def ttest_ind_phenotype(case: Union[SampleCollection, Sample], control: Union[SampleCollection, Sample], alternative="two-sided"):
    """Two sided t-test of case sample(s) and mean expression values in base samples across all genes

        Args:
            case: either Sample of SampleCollection object with case sample(s)

            control: either Sample of SampleCollection object with control sample(s)

            alternative: string with the alternative hypothesis, H1,
                has to be one of the following:
                    - ‘two-sided’: H1: difference in means not equal to value (default)
                    - ‘larger’: H1: difference in means larger than value
                    - ‘smaller’: H1: difference in means smaller than value

        Returns:
            Parameters of t-test
                - tstat: float or numpy array in case of multiple case samples - test statisic
                - pvalue: float or numpy array in case of multiple case samples - pvalue of the t-test
                - df: int or float - degrees of freedom used in the t-test

        """
    l = [np.mean(row) for (idx, row) in control.as_array().iterrows()]
    return ttest_ind(case.as_array(), l, alternative)
