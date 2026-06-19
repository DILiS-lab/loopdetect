import numpy as np
import loopdetect as ld


def test_find_loops_small_matrix():
    jac = np.array([
        [0, 1],
        [1, 0],
    ])

    loops = ld.find_loops(jac)

    assert loops is not None
    assert len(loops) > 0