from patch_analysis import *

def test_area():
    assert(area([[0,0], [0,1], [0,2], [0,3]]) == 4)

def test_perimeter():
    assert(perimeter([[0,0], [0,1], [0,2], [0,3], [1,1]]) == 12)

def test_number_core_areas():
    patch = [[0,0], [0,1], [0,2], [0,3],\
             [1,0], [1,1], [1,2], [1,3],\
             [2,0], [2,1], [2,2], [2,3],\
             [3,1], \
             [4,0], [4,1], [4,2], [4,3],\
             [5,0], [5,1], [5,2], [5,3],\
             [6,0], [6,1], [6,2], [6,3],\
             [7,0], [7,1], [7,2], [7,3]]

    assert(number_core_areas(patch, 1, 60) == 2)
    assert(number_core_areas(patch, 0, 60) == 1)
    assert(number_core_areas(patch, 2, 60) == 0)

def test_get_moore_neighbors():
    moore = get_moore_neighbors([0,0], 60)
    moore.sort()
    expected = [[0,59], [59,0], [1,0], [0,1], [1,59], [59,59], [1,1], [59,1]]
    expected.sort()
    assert(moore==expected)

def test_get_rook_neighbors():
    rook = get_rook_neighbors([0,0], 60)
    rook.sort()
    expected = [[0,59], [59,0], [1,0], [0,1]]
    expected.sort()
    assert(rook==expected)

def test_related_circumscribing_circle():
    pass
