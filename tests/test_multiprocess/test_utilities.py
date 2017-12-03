from multiprocess import available_cores


def test_cores():
    cores_number = available_cores()

    # hope there is always at least one core to use
    assert cores_number >= 1

    # the number of available cores should be integer
    assert type(cores_number) == int
