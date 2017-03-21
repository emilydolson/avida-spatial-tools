from avidaspatial import *


def test_entropy():
    d = {"a": 2, "b": 2}
    assert(np.isclose(entropy(d), 1))
    d = {"a": 3, "b": 6, "c": 1, "d": 2, "e": 1}
    assert(np.isclose(entropy(d), 1.98777337))
    d = {"a": 300, "b": 600, "c": 100, "d": 200, "e": 100}
    assert(np.isclose(entropy(d), 1.98777337))
    d = {"a": 16, "b": 9, "c": 9, "d": 25, "e": 36}
    assert(np.isclose(entropy(d), 2.114357))
    d = {"a": 1, "b": 1, "c": 3, "d": 3, "e": 2, "f": 45}
    assert(np.isclose(entropy(d), 1.0787567))


def test_sqrt_shannon_entropy():
    d = {"a": 1, "b": 1, "c": sqrt(3), "d": sqrt(3), "e": sqrt(2),
         "f": sqrt(45)}
    assert(np.isclose(sqrt_shannon_entropy("tests/grid_task.10000.dat"),
                      entropy(d)))
