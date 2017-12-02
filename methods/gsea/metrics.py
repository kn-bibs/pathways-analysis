from abc import ABC, abstractmethod
from functools import lru_cache

from numpy import mean, std

from utils import jit, abstract_property


# TODO metrics for continuous phenotypes


class DifferentialExpressionMetric(ABC):
    # TODO: metrics could implemented as be standalone functions or builders

    @abstract_property
    def name(self):
        pass

    @abstractmethod
    def __call__(self, expression_profile, control_profile):
        pass


class DifferenceOfClasses(DifferentialExpressionMetric):

    name = 'difference'

    def __call__(self, expression_profile, control_profile):
        return mean(expression_profile) - mean(control_profile)


class RatioOfClasses(DifferentialExpressionMetric):

    name = 'ratio'

    def __call__(self, expression_profile, control_profile):
        return mean(expression_profile) / mean(control_profile)


class SignalToNoise(DifferentialExpressionMetric):

    name = 'signal_to_noise'

    @lru_cache(maxsize=None)
    @jit
    def __call__(self, case, control):
        # TODO check for zero division?

        return (
            (mean(case) - mean(control))
            /
            (std(case) + std(control))
        )


ranking_metrics = {
    metric.name: metric
    for metric in (DifferenceOfClasses, RatioOfClasses, SignalToNoise)
}
