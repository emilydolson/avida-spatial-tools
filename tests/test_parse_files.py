from avidaspatial import *


def test_parse_environment_file_list():
    envs = parse_environment_file_list("tests/example_environment.cfg", (4, 6))
    assert(envs[0].grid == [[set(['test', 'nor']), set(['test', 'nor']),
                             set(['test', 'nor']), set(['test', 'nor'])],
                            [set(['test', 'nor']), set(['nor']), set(['nor']),
                             set([])], [set(['test', 'nor']), set(['nor']),
                            set(['test', 'nor']), set(['test'])],
                            [set(['test', 'nor']), set(['test']), set([]),
                            set(['test'])], [set(['test']), set([]), set([]),
                            set([])], [set([]), set(['test']), set([]), set([])
                                       ]])


def test_load_grid_data():
    data_file = "tests/grid_task.200000.dat"
    data = load_grid_data(data_file)
    assert(data == [[['0b1100'], ['0b0'], ['0b0'], ['0b0'],
                     ['0b0'], ['0b0'], ['0b0'], ['0b0'],
                     ['0b0'], ['0b0'], ['0b0']], [['0b0'],
                    ['0b0'], ['0b0'], ['0b1001'], ['0b1001'],
                    ['0b0'], ['0b0'], ['0b0'], ['0b0'],
                    ['0b0'], ['0b11']], [['0b0'], ['0b0'],
                    ['0b0'], ['0b10'], ['0b1001'], ['0b0'],
                    ['0b0'], ['0b0'], ['0b0'], ['0b0'],
                    ['0b0']], [['0b11'], ['0b0'], ['0b0'],
                    ['0b10'], ['0b0'], ['0b0'], ['0b0'],
                    ['0b0'], ['0b0'], ['0b0'], ['0b0']],
                    [['0b0'], ['0b0'], ['0b0'], ['0b10'],
                     ['0b0'], ['0b0'], ['0b0'], ['0b0'],
                     ['0b0'], ['0b0'], ['0b0']]])


def test_load_grid_data_int():
    data_file = "tests/grid_task.200000.dat"
    data = load_grid_data(data_file, "int")
    assert(data == [[[12], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]],
                    [[0], [0], [0], [9], [9], [0], [0], [0], [0], [0], [3]],
                    [[0], [0], [0], [2], [9], [0], [0], [0], [0], [0], [0]],
                    [[3], [0], [0], [2], [0], [0], [0], [0], [0], [0], [0]],
                    [[0], [0], [0], [2], [0], [0], [0], [0], [0], [0], [0]]])


def test_load_grid_data_list():
    data_files = ["tests/grid_task.300000.dat", "tests/grid_task.200000.dat"]
    data = load_grid_data(data_files)
    assert(data == [[['0b1100', '0b1100'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b0'], ['0b0', '0b0']],

                    [['0b0', '0b0'], ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b1001', '0b1001'], ['0b1001', '0b1001'],
                     ['0b0', '0b0'], ['0b0', '0b0'], ['0b0', '0b1'],
                     ['0b0', '0b0'], ['0b0', '0b0'], ['0b11', '0b11']],

                    [['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b10', '0b10'],
                     ['0b1001', '0b1001'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b1'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b0']],

                    [['0b11', '0b11'], ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b10', '0b10'], ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b1'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b0']],

                    [['0b0', '0b0'], ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b10', '0b10'], ['0b0', '0b0'], ['0b0', '0b0'],
                     ['0b0', '0b0'], ['0b0', '0b1'],
                     ['0b0', '0b0'], ['0b0', '0b0'], ['0b0', '0b0']]])


def test_agg_grid():
    example_grid = [[[1, 2, 1.5], [1, 1, 1]],
                    [[3, 2, 1], [4, 4, 4]]]
    result = [[1.5, 1], [2, 4]]

    agg = agg_grid(example_grid, mean)
    assert(agg == result)
