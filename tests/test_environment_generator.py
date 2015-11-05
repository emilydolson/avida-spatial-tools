import environment_generator as eg
from environment_file_components import *

def test_gen_reaction():
    args = eg.get_args()
    args.maxCount=20
    rxn = gen_reaction(args, "resAND", 1)
    assert(rxn == "REACTION AND and process:resource=resAND:value=2.0:type=pow:frac=0.0025:max=25:depletable=1 requisite:max_count=20\n")

def test_gen_res():
    args = eg.get_args()
    res = gen_res(args, "resAND", 100, .1)
    assert(res == "RESOURCE resAND:inflow=100:outflow=0.1\n")

def test_gen_cell():
    args = eg.get_args()
    args.cellInflow = 10
    args.cellOutflow = 1
    args.inflow = 100
    res = gen_cell(args, "resAND", [1,4,5,6,9,10])
    assert(res == "CELL resAND:1,4,5,6,9,10:inflow=10:outflow=1:initial=100\n")

def test_gen_gradient():
    args = eg.get_args()
    res = gen_gradient(args, "resAND", 100, 5, (4,5))
    assert (res == "GRADIENT_RESOURCE resAND:height=5:plateau=100:spread=4:common=1:updatestep=1000000:peakx=4:peaky=5:plateau_inflow=100:initial=100\n")

if __name__ == "__main__":
    test_gen_reaction()
