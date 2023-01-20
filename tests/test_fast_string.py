from mcdft.fast_string import FastString
from ase.ga.utilities import get_nndist, get_nnmat_string
from mcdft.structure_generator import swap_atoms
import time
from ase.io import read

t1=time.time()
all_configs=[]

atoms=read('POSCAR_AuCuAl')
def get_fast_string_maker(atoms):
    natoms=atoms.copy()
    natoms=natoms*(3,3,3)
    dm = natoms.get_all_distances(mic=True)
    nndist = get_nndist(natoms, dm) + 0.2
    fast_string_maker=FastString(dm,nndist)
    return fast_string_maker


for i in range(100):
    all_configs.append(swap_atoms(atoms.copy()))

fast_string_maker=get_fast_string_maker(atoms)


for config in all_configs:

    #ga_str=get_nnmat_string(config*(3,3,3), mic=True)
    my_str=fast_string_maker.get_nnmat_string(config*(3,3,3))
t2=time.time()
delta_time=t2-t1
print(delta_time)


