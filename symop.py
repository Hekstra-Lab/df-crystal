from os.path import dirname,realpath
import numpy as np
import re

def order_from_rot_mat(mat):
    """Determine the order of a rotation. Returns an integer."""
    det = np.linalg.det(mat)
    tr = int(np.trace(mat))
    if det > 0:
        return [3, 4, 6, 1, 2][tr]
    elif det < 0:
        return [6, 2, 2, 6, 4][tr]
    else:
#TODO: Should this raise an error? Det should be 1 or -1 in theory
        return 1


class symops(dict):
    def __init__(self, libFN=None):
        if libFN is None:
            libFN = dirname(realpath(__file__)) + "/symop.lib"
        self._parse(libFN)
    def _parse(self, libFN):
        with open(libFN, 'r') as f:
            for match in re.findall(r"(?<=\n)[0-9].*?(?=\n[0-9])", '\n'+f.read(), re.DOTALL):
                k = int(match.split()[0])
                self[k] = symop(match)

class symop(dict):
    def __init__(self, text):
        self.number = int(text.split()[0])
        self.name = re.findall(r"(?<=').*?(?=')", text)[0]
        self._parse(text)
    def _parse(self, text):
        for line in text.split('\n')[1:]:
            self[line] = op(line)

class op():
    def __init__(self, text):
        self.rot_mat = np.zeros((3,3))
        ax  = { 
            'X':  np.array([1., 0., 0.]), 
            'Y':  np.array([0., 1., 0.]), 
            'Z':  np.array([0., 0., 1.]),
           }

        for i,t in enumerate(text.split(',')):
            for k,v in ax.items():
                if '-' + k in t:
                    self.rot_mat[:,i] -= v
                elif k in t:
                    self.rot_mat[:,i] += v
        self.rot_mat = self.rot_mat.T #This puts us in accordance with Hovmoller syntax

        self.trans = np.zeros(3)
        div = lambda x: float(x[0])/float(x[1])
        x,y,z = text.split(',')
        self.trans[0] = 0. if '/' not in x else div(re.sub(r"[^\/0-9]", "", x).split('/'))
        self.trans[1] = 0. if '/' not in y else div(re.sub(r"[^\/0-9]", "", y).split('/'))
        self.trans[2] = 0. if '/' not in z else div(re.sub(r"[^\/0-9]", "", z).split('/'))
        self.order = order_from_rot_mat(self.rot_mat)
        o = np.matmul(np.identity(3) - self.rot_mat, np.identity(3) + self.rot_mat)
        self.intrinsic_translation = (1/self.order)*np.matmul(np.sum([np.linalg.matrix_power(self.rot_mat, i) for i in range(self.order)],axis=0), self.trans)

    def transform_xyz(self, vector):
        """
        Transform a 3 vector or nx3 matrix of xyz positions in fractional coordinates.
        """
        return np.matmul(self.rot_mat, vector.T).T + self.trans

    def transform_hkl(self, vector):
        """
        Transform a 3 vector or nx3 matrix of hkl indices.
        """
        return np.matmul(vector, self.rot_mat).astype(int)

    def translate(self, vector):
        """
        There is a decent chance this is garbage. Not necessary now, but TODO: fix this
        """
        return vector + self.trans*vector

class spacegroupnums(dict):
    def __init__(self, libFN=None):
        if libFN is None:
            libFN = dirname(realpath(__file__)) + "/symop.lib"
        self._parse(libFN)
    def _parse(self, libFN):
        with open(libFN, 'r') as f:
            for line in f:
                if line[0] != ' ':
                    k = line.split("'")[1]
                    v = int(line.split()[0])
                    self[k] = v

symops = symops()
spacegroupnums = spacegroupnums()
#These are overrides for CNS files which use a nonstandard nomenclature
spacegroupnums['P 2'] = 3
spacegroupnums['C 2'] = 5
spacegroupnums['P 222'] = 16
spacegroupnums['P 2221'] = 17
spacegroupnums['P 22121'] = 18
spacegroupnums['C 2221'] = 20
spacegroupnums['C 222'] = 21
spacegroupnums['F 222'] = 22
spacegroupnums['I 222'] = 23
spacegroupnums['P 422'] = 89
spacegroupnums['P 421 2'] = 90
spacegroupnums['P 41 22'] = 91
spacegroupnums['P 42 22'] = 93
spacegroupnums['P 43 22'] = 95
spacegroupnums['I 422'] = 97
spacegroupnums['I 41 22'] = 98
spacegroupnums['P 31 21'] = 152
spacegroupnums['P 32 21'] = 154
spacegroupnums['P 622'] = 177
spacegroupnums['P 62 22'] = 180
spacegroupnums['F 41 32'] = 210

spacegroupnames = {v: k for k,v in spacegroupnums.items()}
