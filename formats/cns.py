from io import StringIO
import re
import numpy as np
import pandas as pd
import symop
import crystal

def read(self, cnsfile):
    """
    Initialize attributes and populate the crystal object with data from
    a CNS-formatted reflection file

    Parameters
    ----------
    cnsfile : str or file
        name of a CNS-format file or a file object
    """

    # This could be made fast/lightweight at the expense of readability later
    if isinstance(cnsfile, str):
        cnsfile = open(cnsfile)
    lines = cnsfile.readlines()
    self.header  = [i for i in lines if i[:4] != 'INDE']
    declares     = [i for i in self.header if i[:4] == 'DECL']
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
    usecols  = [1, 2, 3]
    for declare in declares:
        colname = re.search(r"(?<=NAME=)[^\s]*", declare).group()
        if colname != None:
            colnames.append(colname)
            usecols.append(usecols[-1] + 2)
    self.dataname = colnames[-1]


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

def write(self, outfile):
    """
    Write contents of crystal object to a CNS file

    Parameters
    ----------
    outfile : str or file
        name of a CNS file or a file like object
    """
    closeme = False
    if isinstance(outfile, str):
        outfile = open(outfile, 'w')
        closeme = True
        
    outfile.write(''.join(self.header))
    for (h,k,l),d in self.iterrows():
        outfile.write(f"INDEx {h:d} {k:d} {l:d} F= {d['F']:.5f} {d['PHASE']:.7f}\n")

    # If this function opened a file, close it
    if closeme:
        outfile.close()

    return
