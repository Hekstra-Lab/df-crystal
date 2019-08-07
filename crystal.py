from multiprocessing.pool import Pool
from multiprocessing import cpu_count
from io import StringIO
from formats import cns, hkl
import symop, re
import numpy as np
import pandas as pd

def is_phase_key(key):
    """Tell if a string is likely meant to signify a phase in a file header"""
#TODO: This is a bandaid! A systematic technique for identifying phases across input files is essential! Deprecate and replace this ASAP!
    if len(key) < 2:
        return False
    elif key.lower()[:2] == 'ph':
        return True
    else:
        return False


def cellvol(a, b, c, alpha, beta, gamma):
    """
    Compute the volume of a crystallographic unit cell from the lattice constants.

    Parameters
    ----------
    a : float
        Length of the unit cell A axis in angstroms
    b : float
        Length of the unit cell B axis in angstroms
    c : float
        Length of the unit cell C axis in angstroms
    alpha : float
        Unit cell alpha angle in degrees
    beta  : float
        Unit cell beta angle in degrees
    gamma : float
        Unit cell gamma angle in degrees

    Returns
    -------
    float
        The volume of the unit cell in cubic angstroms
    """
    alpha = np.deg2rad(alpha)
    beta  = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)
    V = a*b*c*(1. - np.cos(alpha)*np.cos(alpha) - np.cos(beta)*np.cos(beta) - np.cos(gamma)*np.cos(gamma) + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))**(0.5)
    return V

def orthogonalization(a, b, c, alpha, beta, gamma):
    """
    Compute the orthogonalization matrix from the unit cell constants

    Parameters
    ----------
    a : float
        Length of the unit cell A axis in angstroms
    b : float
        Length of the unit cell B axis in angstroms
    c : float
        Length of the unit cell C axis in angstroms
    alpha : float
        Unit cell alpha angle in degrees
    beta  : float
        Unit cell beta angle in degrees
    gamma : float
        Unit cell gamma angle in degrees

    Returns
    -------
    np.ndarry
        A 3x3 array of the orthogonalization matrix corresponding to the supplied cell parameters
    """
    V = cellvol(a, b, c, alpha, beta, gamma)
    alpha = np.deg2rad(alpha)
    beta  = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)
    O = np.array([
        [a, b*np.cos(gamma), c*np.cos(beta)],
        [0, b*np.sin(gamma), c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)],
        [0., 0., V/(np.sin(gamma)*a*b)]
    ])
    return O

def deorthogonalization(a, b, c, alpha, beta, gamma):
    """
    Compute the deorthogonalization matrix from the unit cell constants

    Parameters
    ----------
    a : float
        Length of the unit cell A axis in angstroms
    b : float
        Length of the unit cell B axis in angstroms
    c : float
        Length of the unit cell C axis in angstroms
    alpha : float
        Unit cell alpha angle in degrees
    beta  : float
        Unit cell beta angle in degrees
    gamma : float
        Unit cell gamma angle in degrees

    Returns
    -------
    np.ndarry
        A 3x3 array of the deorthogonalization matrix corresponding to the supplied cell parameters
    """
    O = orthogonalization(a, b, c, alpha, beta, gamma)
    return np.linalg.inv(O)

def dhkl(h, k, l, a, b, c, alpha, beta, gamma):
    """
    Compute the real space lattice plane spacing, "d", associated with a given hkl.

    Parameters
    ----------
    h : int or np.ndarray
        h-index or indices for which you wish to calculate lattice spacing
    k : int or np.ndarray
        k-index or indices for which you wish to calculate lattice spacing
    l : int or np.ndarray
        l-index or indices for which you wish to calculate lattice spacing
    a : float
        Length of the unit cell A axis in angstroms
    b : float
        Length of the unit cell B axis in angstroms
    c : float
        Length of the unit cell C axis in angstroms
    alpha : float
        Unit cell alpha angle in degrees
    beta  : float
        Unit cell beta angle in degrees
    gamma : float
        Unit cell gamma angle in degrees

    Returns
    -------
    float or array_like
    """
    hkl = np.vstack((h, k, l))
    Oinv = deorthogonalization(a, b, c, alpha, beta, gamma)
    d = 1./np.sqrt(np.sum(np.square(np.matmul(Oinv.T, hkl)), 0))
    return d

class crystalseries(pd.Series):
    spacegroup = None
    cell = None
    A = None
    F = None
    V = None

    @property
    def _constructor(self):
        return crystalseries

class crystal(pd.DataFrame):
    """
    Representation of a crystal

    Attributes
    ----------
    spacegroup : int
        Number corresponding to the crystal space group
    cell : np.ndarray
        Unit cell constants of crystal (A, B, C, alpha, beta, gamma)
    A : np.ndarray
        Matrix of unit cell vectors
    V : float
        Unit cell volumen in cubic angstroms
    """

    _metadata = ['header', 'spacegroup', 'cell', 'A', 'V']
    header = None
    spacegroup = None
    cell = None
    A = None
    V = None
    datacol = None
    errorcol= None

    @property
    def _constructor(self):
        return crystal

    @property
    def _constructor_sliced(self):
        return crystalseries

    def read_cns(self, cnsfile):
        """
        Initialize attributes and populate the crystal object with data
        from a cns formatted reflection file
        
        Parameters
        ----------
        cnsfile : str or file
            name of an hkl file or a file like object
        """
        return cns.read(self, cnsfile)

    def write_cns(self, outfile):
        """
        Write contents of crystal object to a CNS file

        Parameters
        ----------
        outfile : str or file
            name of an hkl file or a file like object
        """
        return cns.write(self, outfile)

    def read_hkl(self, hklfile, a=None, b=None, c=None, alpha=None,
                 beta=None, gamma=None, sg=None):
        """
        Initialize attributes and populate the crystal object with data 
        from a HKL file of reflections. This is the output format used 
        by Precognition when processing Laue diffraction data.

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
        return hkl.read(self, hklfile, a, b, c, alpha, beta, gamma, sg)
        
    def _label_centrics(self):
        """
        Add 'CENTRIC' key to self. Label centric reflections as True.
        """
        self['CENTRIC'] = False
        hkl = np.vstack(self.index)
        centric = np.zeros(len(hkl), dtype=bool)
        for key,op in symop.symops[self.spacegroup].items():
            centric = np.all(np.isclose(op(hkl.T), -hkl.T), 0) | centric
        self['CENTRIC'] = centric
        return self

    def remove_centrics(self):
        """
        Return a copy of self without Miller indices of centric reflections.

        Returns
        -------
        crystal : crystal
            Copy of self without centric reflections
        """
        hkl = np.vstack(self.index)
        centric = np.zeros(len(hkl), dtype=bool)
        for key,op in symop.symops[self.spacegroup].items():
            centric = np.all(np.isclose(op(hkl.T), -hkl.T), 0) | centric
        return self[~centric]

    def hkl_to_reciprocal_asu(self, hkl):
        hkl = np.array(hkl)
        labels = None
        for key,op in symop.symops[self.spacegroup].items():
            if labels is None:
                labels = op(hkl.T)[None, :, :]
            else:
                labels = np.concatenate((labels, op(hkl.T)[None, :, :]), 0)
        labels = np.sort(labels, 0)[-1].astype(int)
        merged = crystal(index=self.index)
        print(labels.shape)
        merged['MERGEDH'], merged['MERGEDK'], merged['MERGEDL'] = labels
        return merged

    def populate_merged_hkls(self):
        hkl = np.vstack(self.index)
        merged = self.hkl_to_reciprocal_asu(hkl)
        #self.update(merged)
        self['MERGEDH'] = merged['MERGEDH']
        self['MERGEDK'] = merged['MERGEDK']
        self['MERGEDL'] = merged['MERGEDL']
        self._coerce_dtypes()
        return self

    def _coerce_dtypes(self):
        """
        This needs to exist to correct some problematic default behaviors in pandas. 
        In future releases, this will hopefull be unnecessary. 
        """
        indexnames = self.index.names
        self.reset_index(inplace=True)
        datatypes = {
            'MERGEDH' : int , 
            'MERGEDK' : int , 
            'MERGEDL' : int , 
            'H' : int , 
            'K' : int , 
            'L' : int , 
        }
        for k,v in datatypes.items():
            if k in self:
                self[k] = self[k].astype(v)
        self.set_index(indexnames, inplace=True)
        return self

    def unmerge_anomalous(self):
        F = self.reset_index()
        F['MERGEDH'] = F['H']
        F['MERGEDK'] = F['K']
        F['MERGEDL'] = F['L']
        Friedel = F.copy()
        Friedel[['H', 'K', 'L']] = -Friedel[['H', 'K', 'L']]
        for k in Friedel:
            if is_phase_key(k):
                Friedel[k] = np.rad2deg(np.angle(np.exp(1j*np.pi)*np.exp(1j*np.deg2rad(Friedel[k]))))
        F = F.append(Friedel).set_index(['H', 'K', 'L'])
        F = F[~F.index.duplicated(keep='first')]
        #self.update(F)
        self._coerce_dtypes()
        self.__init__(F)
        return self

    def unmerge(self):
        F = self.reset_index()
        F['MERGEDH'] = F['H']
        F['MERGEDK'] = F['K']
        F['MERGEDL'] = F['L']
        for k,op in symop.symops[self.spacegroup].items():
            f = self.reset_index().copy()
            f['MERGEDH'] = f['H']
            f['MERGEDK'] = f['K']
            f['MERGEDL'] = f['L']
            f[['H', 'K', 'L']] = np.array(op(f[['H', 'K', 'L']].T).T, int)
            F = F.append(f)
        Friedel = F.copy()
        Friedel[['H', 'K', 'L']] = -Friedel[['H', 'K', 'L']]
        for k in Friedel:
            if is_phase_key(k):
                Friedel[k] = np.rad2deg(np.angle(np.exp(1j*np.pi)*np.exp(1j*np.deg2rad(Friedel[k]))))
        F = F.append(Friedel).set_index(['H', 'K', 'L'])
        F = F[~F.index.duplicated(keep='first')]
        self._coerce_dtypes()
        self.__init__(F)
        return self

    def rotate(self, phistep, axis=None):
        """
        Rotate the crystal unit cell by phistep about an axis. Update the A matrix of the crystal accordingly. 

        Parameters
        ----------
        phistep : float
            The number of degrees to rotate the crystal
        axis : np.ndarray
            The cartesian axis on which to rotate the crystal

        Returns
        -------
        crystal
            This method returns self for easy chaining. 
        """
        phistep = np.deg2rad(phistep)
        if axis is None:
            axis = np.array([0., 1, 0.])
        #Formula for arbitrary rotation about an axis
        R = np.identity(3)*np.cos(phistep) + np.sin(phistep)*np.cross(axis, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]).T + (1 - np.cos(phistep))*np.outer(axis, axis)
        self.A = np.matmul(self.A, R)
        return self

    def reflections(self, wavelength=None, tol=None, detector_distance=None):
        """
        Parameters
        ----------
        wavelength : float
            The wavelength of the x-ray beam in Angstroms
        tol : float, optional
            Allowed error in the Bragg condition to accept a reflection. The defalt value is 0.001.
        detector_distance : float
            The distance of the simulated "detector" from the crystal position. The default value is 100 mm.

        Returns
        -------
        pd.Dataframe
            Dataframe object with reflections, structure factors. 
        """
        detector_distance= 100. if detector_distance is None else detector_distance
        wavelength = 1. if wavelength is None else wavelength
        tol = 0.001 if tol is None else tol
        Ainv = np.linalg.inv(self.A)
        So = np.array([0, 0, 1./wavelength])
        def err(x):
            h = np.array(x.name)
            S = np.matmul(Ainv, h)
            return 0.5*np.dot(S, S) - np.dot(S, So)
        F = self[np.abs(self.apply(err, 1)) <= tol]
        def coordinate(x):
            h = np.array(x.name)
            S = np.matmul(Ainv, h)
            S1 = S+So
            S1 = S1/(wavelength*np.linalg.norm(S1)) #Map to Ewald Sphere
            XYZ = detector_distance*S1/S1[2]
            return pd.Series(XYZ)
        F['A'] = [self.A[:,0]]*len(F)
        F['B'] = [self.A[:,1]]*len(F)
        F['C'] = [self.A[:,2]]*len(F)
        return F.join(F.apply(coordinate, 1).rename(columns={0:'X', 1:'Y', 2:'Z'}))

    def orientrandom(self):
        """
        Randomly rotate the unit cell and update the A-matrix correspondingly. 
#TODO: this implementation doesn't give uniform sampling! reimplement this using quaternions as in https://marc-b-reynolds.github.io/distribution/2017/01/27/UniformRot.html

        Returns
        -------
        crystal
            This method resturns self for easy chaining. 
        """
        self.rotate(360.*np.random.random(), axis=[1., 0., 0.])
        self.rotate(360.*np.random.random(), axis=[0., 1., 0.])
        self.rotate(360.*np.random.random(), axis=[0., 0., 1.])
        return self

    def phiseries(self, phistep, nsteps, reflections_kwargs=None, axis=None, nprocs=None):
        """
        Compute a series of images by rotating the crystal. This method uses multiprocessing for parallelization. 

        Parameters
        ----------
        phistep : float
            Phi angle step in degrees between frames
        nsteps : int
            Number of images to simulate
        reflections_kwargs : dict
            Keword arguments to pass to crystal.reflections in case you want to override the defaults. Default is None.
        axis : np.ndarray
            Axis about which to rotate the crystal. 
        nprocs : int
            Number of processors to use for this calculation. 

        Returns
        -------
        pd.DataFrame
            Datframe containing the accepted reflections from the rotation series. 
        """
        axis = [0,1,0] if axis is None else axis
        reflections_kwargs = {} if reflections_kwargs is None else reflections_kwargs
        iterable = [(self.copy().rotate(i*phistep, axis=axis), reflections_kwargs, i) for i in range(nsteps)]
        #iterable = [i*phistep for i in range(nsteps)]
        nprocs = cpu_count() if nprocs is None else nprocs
        p = Pool(nprocs)
        F = p.map(_phihelper, iterable)
        p.terminate()
        return pd.concat(F)

def _phihelper(X):
    cryst, kw, i = X
    f = cryst.reflections(**kw)
    f['PHINUMBER'] = i
    return f

