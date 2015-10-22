from parse_files import *

def test_parse_environment_file():
    env = parse_environment_file("example_environment.cfg", (4,6))
    assert(env == [[set(['test', 'nor']), set(['test', 'nor']), \
                    set(['test', 'nor']), set(['test', 'nor'])],\
                   [set(['test', 'nor']), set(['nor']), set(['nor']), \
                    set([])], [set(['test', 'nor']), set(['nor']), \
                    set(['test', 'nor']), set(['test'])],\
                   [set(['test', 'nor']), set(['test']), set([]), \
                    set(['test'])], [set(['test']), set([]), set([]), \
                    set([])], [set([]), set(['test']), set([]), set([])]])

def test_load_grid_data():
    data_file = "grid_task.200000.dat"
    data = load_grid_data(data_file, 5)
    assert(data == [[['0b00000'], ['0b00000'], ['0b00000'], ['0b00000'],\
                     ['0b00000'], ['0b00000'], ['0b00000'], ['0b00000'], \
                     ['0b00000'], ['0b00000'], ['0b00000']], [['0b00000'],\
                    ['0b00000'], ['0b00000'], ['0b01001'], ['0b01001'],\
                    ['0b00000'], ['0b00000'], ['0b00000'], ['0b00000'],\
                    ['0b00000'], ['0b00011']], [['0b00000'], ['0b00000'],\
                    ['0b00000'], ['0b00010'], ['0b01001'], ['0b00000'], \
                    ['0b00000'], ['0b00000'], ['0b00000'], ['0b00000'], \
                    ['0b00000']], [['0b00011'], ['0b00000'], ['0b00000'],\
                    ['0b00010'], ['0b00000'], ['0b00000'], ['0b00000'],\
                    ['0b00000'], ['0b00000'], ['0b00000'], ['0b00000']],\
                    [['0b00000'], ['0b00000'], ['0b00000'], ['0b00010'],\
                     ['0b00000'], ['0b00000'], ['0b00000'], ['0b00000'],\
                     ['0b00000'], ['0b00000'], ['0b00000']]])

def test_agg_grid():
    example_grid = [ [ [1,2,1.5], [1,1,1] ],\
                     [ [3,2,1], [4,4,4] ] ]
    result = [[1.5, 1], [2, 4]]

    agg = agg_grid(example_grid, mean)
    assert(agg == result)
