from numpy.ma import sqrt
from pytest import approx

from metrics import difference_of_classes, signal_to_noise, ratio_of_classes


def test_difference():
    assert difference_of_classes([10], [5]) == 5
    assert difference_of_classes([1, 2, 3], [0, 1, 2]) == 1


def test_signal_to_noise():
    assert signal_to_noise((1, 2, 3), (0, 1, 2)) == approx((2 - 1) / (sqrt(1) + sqrt(1)))


def test_ratio_of_classes():
    assert ratio_of_classes([1, 2, 3], [0, 1, 2]) == 2