from io import StringIO
import re
import numpy as np
import pandas as pd
import symop
import crystal

def read(self, hklfile, a=None, b=None, c=None, alpha=None, beta=None,
         gamma=None, sg=None):
    """
    Initialize attributes and populate the crystal object with data from
    a HKL file of reflections. This is the output format used by 
    Precognition when processing Laue diffraction data.

    Parameters
    ----------
    hklfile : str or file
        name of an hkl file or a file object
    a : float
        edge length, a, of the unit cell
    b : float
        edge length, b, of the unit cell
    c : float
        edge length, c, of the unit cell
    alpha : float
        interaxial angle, alpha, of the unit cell
    beta : float
        interaxial angle, beta, of the unit cell
    gamma : float
        interaxial angle, gamma, of the unit cell
    sg : str or int
        If int, this should specify the space group number. If str, 
        this should be a space group symbol
    """
    # Read data from HKL file
    usecols = [0, 1, 2, 3, 4]
    F = pd.read_csv(hklfile, header=None, delim_whitespace=True,
                       names=["H", "K", "L", "F", "SigF"],
                       usecols=usecols)
    for k,v in F.items():
        self[k] = v

    self.set_index(["H", "K", "L"], inplace=True)

    # Set Crystal attributes
    if (a and b and c and alpha and beta and gamma):
        self.cell = np.array([a, b, c, alpha, beta, gamma])
        self.A = crystal.orthogonalization(*self.cell).T
        self.V = crystal.cellvol(*self.cell)
    if sg:
        sg = re.sub(r'\(', '', sg)
        sg = re.sub(r'\)', ' ', sg)
        sg = sg[0] + ' ' + sg[1:].strip()
        self.spacegroup = symop.spacegroupnums[sg]
        
    return self

def write(self, outfile, sf_key="F", err_key="SigF", phase_key=None,
          weight_key=None):
    """
    Write contents of crystal object to an HKL file

    Parameters
    ----------
    outfile : str or file
        name of an hkl file or file-like object
    sf_key : str
        key for structure factor in DataFrame
    err_key : str
        key for structure factor error in DataFrame
    phase_key : str
        key for phase in DataFrame
    weight_key : str
        key for structure factor weights in DataFrame
    """
    closeme = False
    if isinstance(outfile, str):
        outfile = open(outfile, 'w')
        closeme = True
        
    for (h,k,l), d in self.iterrows():
        if phase_key is None and weight_key is None:
            outfile.write("{h:5d}{k:5d}{l:5d}{d[sf_key]:15.2f}{d[err_key]:15.2f}\n")
        elif phase_key and weight_key is None:
            outfile.write("{h:5d}{k:5d}{l:5d}{d[sf_key]:15.2f}{d[err_key]:15.2f}{d[phase_key]:15.7f}\n")
        else:
            outfile.write("{h:5d}{k:5d}{l:5d}{d[sf_key]:15.2f}{d[err_key]:15.2f}{d[phase_key]:15.7f}{d[weight_key]:15.7f}\n")

    # If this function opened a file, close it
    if closeme:
        outfile.close()
            
    return
