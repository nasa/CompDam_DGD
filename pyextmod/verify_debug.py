import numpy as np
import CompDam_DGD
import helpers as h
import os, sys, argparse, inspect, json, re, numpy, shutil

sv_attributes_ignore = ('old', 'debugpy_count', 'direct')

def _main_entry(args):
    '''
    Loads debug file
    Runs compdam using sv_old
    Writes a json file with the state variables calculated by Abaqus and the Python Extension Module
    '''

    # debug file name
    debug_file_name = os.path.basename(args.debugfile)
    debug_file_name_no_ext = os.path.splitext(debug_file_name)[0]
    # location to path
    debug_file_abspath = os.path.abspath(args.debugfile)
    sys.path.append(os.path.dirname(debug_file_abspath))
    # job name
    jobName = re.sub(r'-[0-9]*-debug-[0-9]*', '', debug_file_name_no_ext)

    # For debugging
    with open(os.path.join(os.path.dirname(debug_file_abspath), jobName+'_python.log'), 'w') as log:
        log.write('Debug file name: ' + debug_file_name + '\n')

        # logging
        pyextmod_log_file_name = jobName+'_fortran.log'
        CompDam_DGD.dgd_mod.log_init(level=4, filename=pyextmod_log_file_name, totaltime=0.081)

        # Load the debug file
        (m, p, sv, sv_old, F, F_old, U, debugpy) = h.loaddebugpy(filename=debug_file_name_no_ext)
        log.write('\nState variables: \n' + str(sv) + '\n')
        print(sv)

        # Run CompDam
        sv_calculated  = sv_old
        Cauchy = np.zeros((3,3), order='F')
        enerintern = 0
        enerinelas = 0
        func = getattr(CompDam_DGD.dgd_mod, args.subroutine)
        func(u=U, f=F, f_old=F_old, m=m, p=p, sv=sv_calculated, ndir=3, nshr=3, dt=0, density_abq=debugpy.density_abq, cauchy=Cauchy, enerintern=enerintern, enerinelas=enerinelas)

        # Move the pyextmod log file to the testoutput directory
        CompDam_DGD.dgd_mod.log_close()
        os.rename(os.path.abspath(pyextmod_log_file_name), os.path.abspath(os.path.join(os.pardir, 'tests', 'testOutput', pyextmod_log_file_name)))

        ## GOAL is to compare state variables --> if state variables are computed correctly, assume the debug.py file logic works
        sv_comparison = {}
        attributes = _get_attributes(sv)
        for a in attributes:
            if a in sv_attributes_ignore:
                continue
            abq_sv = getattr(sv, a)
            pyextmod_sv = getattr(sv_calculated, a)
            # Convert numpy arrays to lists for serialization
            if isinstance(abq_sv, numpy.ndarray):
                abq_sv = abq_sv.tolist()
            if isinstance(pyextmod_sv, numpy.ndarray):
                pyextmod_sv = pyextmod_sv.tolist()
            sv_comparison[a] = (abq_sv, pyextmod_sv)

        # Write to file
        output_filename = jobName +'_pyextmod_results.json'
        with open(os.path.join(os.path.dirname(debug_file_abspath), output_filename), 'w') as outfile:
            json.dump(sv_comparison, outfile, indent=2)

        log.write('End of python extension module execution\n')

    return


def _get_attributes(obj):
    attributes = []
    for i in inspect.getmembers(obj):
        if not i[0].startswith('_'):
            if not inspect.ismethod(i[1]):
                attributes.append(i[0])
    return attributes


if __name__ == "__main__":
	# Arguments
    parser = argparse.ArgumentParser(description='Loads in state variables from a debug file to prove continuity between debug file and CompDam code.')
    parser.add_argument('subroutine', choices=['dgdinit', 'dgdevolve', 'dgdkinkband'], help='Specify which functionality from DGD to call.')
    parser.add_argument('debugfile', action='store', help='Path to debug file.')
    parser.set_defaults(func=_main_entry)
    
    # Parse the args
    args = parser.parse_args()
    args.func(args)
