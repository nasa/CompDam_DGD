import numpy as np
import CompDam_DGD
import os, sys


def loaddebugpy(filename='', jobname='', suffix=''):
    '''
    This function provides functionality to load a debug.py file
    Specifiy filename or jobname
    Assumes debug.py is in working directory
    '''

    # Handle arguments
    if (filename and jobname) or (not filename and not jobname):
        raise Exception('Must specify either filename or jobname')
    if jobname:
        filename = jobname+'-debug'

    # Import the debug file
    debugpy = __import__(filename+suffix)

    # Material properties
    first_line = [int(debugpy.featureFlags['integer']), 2., debugpy.thickness, 0., 0., 0., 0., 0., ]
    m = CompDam_DGD.matprop_mod.loadmatprops('IM7-8552', 40, first_line + debugpy.m)
    CompDam_DGD.matprop_mod.consistencychecks(m, issuewarnings=False)

    # Parameters
    p = debugpy.p
    default_parameters = CompDam_DGD.parameters_mod.loadparameters()
    for k, v in p.items():
        setattr(default_parameters, k.lower(), v)
    p = default_parameters
    print(p)

    # Deformations
    F = tensorAsListTo3x3(debugpy.F)
    F_old = tensorAsListTo3x3(debugpy.F_old)
    U = tensorAsListTo3x3(debugpy.U)

    # State variables
    sv = CompDam_DGD.statevar_mod.loadstatevars(len(debugpy.sv), debugpy.sv, m)
    sv_old = CompDam_DGD.statevar_mod.loadstatevars(len(debugpy.sv_old), debugpy.sv_old, m)

    return (m, p, sv, sv_old, F, F_old, U)


def tensorAsListTo3x3(t):
    '''
    Convert a list (len=9) to a 3x3 numpy array
    Follows abaqus VUMAT component ordering
    '''
    n = np.array([
        [t[0], t[3], t[8]],
        [t[6], t[1], t[4]],
        [t[5], t[7], t[2]] ])
    return n