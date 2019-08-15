from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import crystal
from subprocess import call
from os import devnull
from io import StringIO

pdbs = pd.read_csv(StringIO("""PDBID, space group number, space group name
6gl4, 3, P 1 2 1
6ofl, 4, P 1 21 1        
6h7c, 5, C 1 2 1
2z51, 16, P 2 2 2
6a8k, 17, P 2 2 21
6nsv, 18, P 2 21 21
3kxe, 19, P 21 21 21     
1oel, 20, C 2 2 21       
6E1X, 21, C 2 2 2
6fxw, 22, F 2 2 2
3ruw, 23, I 2 2 2        
5CR4, 24, I 21 21 21
2ON8, 75, P 4
6E6T, 76, P 41
4AK8, 77, P 42
6E6N, 78, P 43
3t4d, 79, I 4            
6GUS, 80, I 41
5ZOA, 89, P 4 2 2
5Z3A, 90, P 4 21 2
5QRH, 91, P 41 2 2
6H9J, 92, P 41 21 2
6Q8D, 93, P 42 2 2
6O11, 94, P 42 21 2
6OH9, 95, P 43 2 2
9lyz, 96, P 43 21 2      
6NMT, 97, I 4 2 2
6NRH, 98, I 41 2 2
6mbu, 143, P 3           
6b8z, 152, P 31 2 1      
5w79, 154, P 31 1 2      
6e02, 168, P 6           
6ovt, 169, P 61          
6cy6, 177, P 6 2 2       
6h9p, 180, P 62 2 2      
4cy9, 210, F 41 3 2
"""
), dtype={"PDBID": str, "space group number": int, "space group name": str})



def test(pdbid, high_resolution=4., verbose=False):
    null = open(devnull, 'w')
    call(f'wget files.rcsb.org/download/{pdbid}.pdb'.split(), stdout=null, stderr=null)
    call(f'phenix.fmodel {pdbid}.pdb high_resolution={high_resolution}'.split(), stdout=null)
    call(f'phenix.reflection_file_converter {pdbid}.pdb.mtz --cns={pdbid}.cns'.split(), stdout=null)
    call(f'phenix.reflection_file_converter {pdbid}.pdb.mtz --expand_to_p1 --cns={pdbid}_p1.cns'.split(), stdout=null)
    null.close()

    #P1 control data file
    p1_filename = f"{pdbid}_p1.cns" 
    merged_filename = f"{pdbid}.cns" 

    ctrl = crystal.Crystal().read_cns(p1_filename)

    test_in = crystal.Crystal().read_cns(merged_filename)
    test = test_in.copy().unmerge()

    passed = True

    if verbose:
        print("Checking structure factor amplitudes...", end=' ')
    check = np.all(np.isclose(test.F, ctrl.loc[test.index, 'F']))
    passed = check & passed
    if verbose:
        if np.all(check):
            print("OK")
        else:
            print("Fail")

    if verbose:
        print("Plotting reciprocal space coverage...")
        test.plot_reciprocal_coverage(bins=200)

    if verbose:
        print("Checking phases...", end=' ')
    problematic = test.loc[~np.isclose(
        np.exp(1j*np.deg2rad(ctrl.loc[test.index].PHASE)), 
        np.exp(1j*np.deg2rad(test.PHASE)), 
        )
    ]
    good = test.loc[np.isclose(
        np.exp(1j*np.deg2rad(ctrl.loc[test.index].PHASE)), 
        np.exp(1j*np.deg2rad(test.PHASE)), 
        )
    ]

    if len(problematic) == 0:
        if verbose:
            print("OK")
    else:
        passed = False
        if verbose:
            print("Fail")
            print(f"{100*len(problematic)/len(test):0.2f} %% incorrect phases")
            centric = problematic[problematic.CENTRIC]
            acentric = problematic[~problematic.CENTRIC]
            print(f"{100*problematic.CENTRIC.sum()/len(problematic):0.2f} %% centrics")
            print(pd.crosstab(centric.op, centric.Friedel))
            print(f"Plotting phase error for centric reflections... ")
            plt.hexbin(ctrl.loc[centric.index, 'PHASE'], centric.PHASE)
            plt.show()

            print(f"Plotting phase error for acentric reflections... ")
            print(pd.crosstab(acentric.op, acentric.Friedel))
            plt.hexbin(ctrl.loc[acentric.index, 'PHASE'], acentric.PHASE)
            plt.show()
            acentric.plot_reciprocal_coverage(bins=200)
            acentric.plot_reciprocal_coverage_3d()

    call(f'rm {pdbid}.pdb {pdbid}.cns {pdbid}_p1.cns {pdbid}.pdb.mtz'.split())
    return passed

if __name__=="__main__":
    from sys import argv
    if len(argv) == 1:
        pdbs['Test'] = pdbs['PDBID'].apply(test)
        print(pdbs)
    else:
        for pdbid in argv[1:]:
            test(pdbid, verbose=True)

