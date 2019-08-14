from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import crystal
from subprocess import call
from os import devnull


def test(pdbid, high_resolution=4.):
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


    print("Checking structure factor amplitudes...", end=' ')
    check = np.isclose(test.F, ctrl.loc[test.index, 'F'])
    if np.all(check):
        print("OK")
    else:
        print("Fail")

    print("Plotting reciprocal space coverage...")
    test.plot_reciprocal_coverage(bins=200)

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
        print("OK")
    else:
        print("Fail")
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

if __name__=="__main__":
    from sys import argv
    test(argv[1])

