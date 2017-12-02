from multiprocess import Pool


def test_imap():
    data = [1, 2, 3, 4]

    # for single process Pool has simplified logic,
    #  thus requires separate test iteration
    for num_processes in [1, 2]:

        pool = Pool(num_processes)

        # let's square: func = pow(x, 2)
        squared = pool.imap(pow, data, shared_args=[2])

        assert set(squared) == {1, 4, 9, 16}
