import os
import tempfile
import subprocess

def write(self, outfile, sf_key="F", err_key="SigF", phase_key=None,
          weight_key=None):
    """
    Write contents of crystal object to an MTZ file

    Parameters
    ----------
    outfile : str
        filename of an MTZ file to write
    sf_key : str
        key for structure factor in DataFrame
    err_key : str
        key for structure factor error in DataFrame
    phase_key : str
        key for phase in DataFrame
    weight_key : str
        key for structure factor weights in DataFrame
    """
    # Write out temporary HKL file
    tmpfp = tempfile.NamedTemporaryFile("w", suffix=".hkl", delete=False)
    self.write_hkl(tmpfp, sf_key, err_key, phase_key, weight_key)
    tmpfp.close()
    
    # Build up call to f2mtz
    command = f"f2mtz HKLIN {tmpfp.name} HKLOUT {outfile} << EOF\n"
    a, b, c, alpha, beta, gamma = tuple(self.cell)
    command += f"CELL {a} {b} {c} {alpha} {beta} {gamma}\n"
    
    if phase_key is None and weight_key is None:
        command += "LABOUT H   K  L   FP SIGFP\n"
        command += "CTYPE  H   H  H   F  Q\n"
        command += "FORMAT '(3(F5.0),2(F15.3))'\n"
    elif phase_key and weight_key is None:
        command += "LABOUT H   K  L   FP SIGFP PHIC\n"
        command += "CTYPE  H   H  H   F  Q     P\n"
        command += "FORMAT '(3(F5.0),2(F15.3),F15.7)'\n"
    else:
        command += "LABOUT H   K  L   FP SIGFP PHIC WT\n"
        command += "CTYPE  H   H  H   F  Q     P    W\n"
        command += "FORMAT '(3(F5.0),2(F15.3),2(F15.7))'\n"
        
    command += f"SYMM {self.spacegroup}\n"
    command += "EOF"

    # Call f2mtz
    p = subprocess.call(command, stdout=subprocess.DEVNULL, shell=True)

    # Remove temporary file
    os.remove(tmpfp.name)

    return
