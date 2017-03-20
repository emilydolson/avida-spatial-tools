from avidaspatial import *


def test_mode():
    assert(mode([1, 1, 1, 2, 3]) == 1)
    assert(mode([5, 10, 3, 2, 5, 1, 11]) == 5)


def test_mean():
    assert(mean([1, 2, 3]) == 2)
    assert(mean([2, 2]) == 2)
    assert(mean([5, 3, 10]) == 6)
    assert(mean([1, 2]) == 1.5)


def test_median():
    assert(median([1, 5, 3]) == 3)


def test_string_avg():
    strings = ["0b0", "0b10", "0b011"]
    assert(string_avg(strings) == "0b010")


def test_convert_to_pysal():
    data = [[1, 2], [3, 4]]
    w, data = convert_to_pysal(data)
    assert(w[0] == {1: 1.0, 2: 1.0})
    assert(data[0] == [1])


def test_get_world_dimensions():
    dims = get_world_dimensions("grid_task.100000.dat")
    assert(dims == (11, 5))
