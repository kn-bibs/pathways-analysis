from abc import ABC, abstractmethod
from functools import lru_cache

from numpy import mean, std

from utils import jit

# TODO metrics for continuous phenotypes


class DifferentialExpressionMetric(ABC):
    # TODO: metrics could implemented as be standalone functions or builders

    @abstractmethod
    def __call__(self, expression_profile, control_profile):
        pass


class DifferenceOfClasses(DifferentialExpressionMetric):

    def __call__(self, expression_profile, control_profile):
        return mean(expression_profile) - mean(control_profile)


class RatioOfClasses(DifferentialExpressionMetric):

    def __call__(self, expression_profile, control_profile):
        return mean(expression_profile) / mean(control_profile)


class SignalToNoise(DifferentialExpressionMetric):

    @lru_cache(maxsize=None)
    @jit
    def __call__(self, case, control):
        # TODO check for zero division?

        return (
            (mean(case) - mean(control))
            /
            (std(case) + std(control))
        )

