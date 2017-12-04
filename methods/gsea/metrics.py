# TODO metrics for continuous phenotypes
from inspect import signature
from functools import lru_cache
from typing import Iterable

from numpy import mean
from statistics import stdev
# for population/biased standard deviation use:
# from numpy import std as stdev

from utils import jit


RANKING_METRICS = {}


def differential_expression_metric(case: Iterable[float], control: Iterable[float]):
    """Template for differential-expression metrics.

    Args:
        case: expression profile of the first class of samples
        control:  expression profile of the second class of samples
    """
    pass


def metric(name):
    """Decorates differential-expression metric.

    Args:
        name: user-visible name of the metric
    """
    def decorator(func):
        func.name = name

        func_signature = signature(func)
        typed_signature = signature(differential_expression_metric)

        # check only names and count of parameters
        if func_signature.parameters.keys() != typed_signature.parameters.keys():
            raise NameError(
                f'Signature of "{name}" metric does not match '
                f'the template: {typed_signature}'
            )

        # replace signature to have type annotation
        func.__signature__ = typed_signature

        # save the metric
        RANKING_METRICS[name] = func
        return func

    return decorator


@metric('difference')
def difference_of_classes(case, control):
    return mean(case) - mean(control)


@metric('ratio')
def ratio_of_classes(case, control):
    return mean(case) / mean(control)


@metric('signal_to_noise')
@lru_cache(maxsize=None)
@jit
def signal_to_noise(case, control):
    """Calculates SNR as ratio of means difference and deviation sum.

    Case and control has to be tuples or other hashable iterable.

    Assumes that there are:
        - at least two samples in both case and control
        - the samples have non-zero variation
    """
    return (
        (mean(case) - mean(control))
        /
        (stdev(case) + stdev(control))
    )
