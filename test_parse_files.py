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

