from io import StringIO
import re
import numpy as np
import pandas as pd
import symop
import crystal

def read_cns(self, hklfile):
    """
    Initialize attributes and populate the crystal object with data from
    a cns formatted reflection file

    Parameters
    ----------
    hklfile : str or file
        name of an hkl file or a file like object
    """

    # This could be made fast/lightweight at the expense of readability later
    if isinstance(hklfile, str):
        hklfile = open(hklfile)
    lines = hklfile.readlines()
    self.header  = [i for i in lines if i[:4] != 'INDE']
    declare      = [i for i in self.header if i[:4] == 'DECL'][0]
    lines        = [i for i in lines if i[:4] == 'INDE']

    a = float(re.search(r'(?<=a=)[^\s]+(?<!\s)', ''.join(self.header)).group())
    b = float(re.search(r'(?<=b=)[^\s]+(?<!\s)', ''.join(self.header)).group())
    c = float(re.search(r'(?<=c=)[^\s]+(?<!\s)', ''.join(self.header)).group())
    alpha = float(re.search(r'(?<=alpha=)[^\s]+(?<!\s)', ''.join(self.header)).group())
    beta  = float(re.search(r'(?<=beta=)[^\s]+(?<!\s)', ''.join(self.header)).group())
    gamma = float(re.search(r'(?<=gamma=)[^\s]+(?<!\s)', ''.join(self.header)).group())
    self.cell = np.array([a, b, c, alpha, beta, gamma])
    self.A = crystal.orthogonalization(*self.cell).T
    self.V = crystal.cellvol(*self.cell)

    sg = re.search(r'(?<=sg=)[^\s]+(?<!\s)', ''.join(self.header)).group()
    if sg == "P2(1)":
        self.spacegroup = 4
    else:
        sg = re.sub(r'\(', '', sg)
        sg = re.sub(r'\)', ' ', sg)
        sg = sg[0] + ' ' + sg[1:].strip()
        self.spacegroup = symop.spacegroupnums[sg]

    colnames = ['H', 'K', 'L']
    colnames.append(re.search(r"(?<=NAME=)[^\s]*", declare).group())
    self.dataname = colnames[-1]

    usecols  = [1, 2, 3, 5]

    # Determine if there is phase information in the file
    if len(lines[0].split()) == 7:
        colnames.append('PHASE')
        usecols.append(6)

    f = StringIO(''.join(lines))
    F = pd.read_csv(f, 
                    delim_whitespace=True, 
                    names=colnames,
                    usecols=usecols,
    )

    for k,v in F.items():
        self[k] = v

    self['D'] = crystal.dhkl(F['H'], F['K'], F['L'], *self.cell)
    self.set_index(['H', 'K', 'L'], inplace=True)
    self._label_centrics()
    return self

def write_cns(self, outfile):
    """
    Write contents of crystal object to a CNS file

    Parameters
    ----------
    outfile : str or file
        name of an hkl file or a file like object
    """
    if isinstance(outfile, str):
        outfile = open(outfile, 'w')

    outfile.write(''.join(self.header))
    for (h,k,l),d in self.iterrows():
        outfile.write("".format(h, k, l, self.datacol))
