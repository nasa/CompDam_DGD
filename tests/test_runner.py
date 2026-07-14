#
# Unittest code to run tests on single element models with DGD
#

import os
import shutil
import sys
import abaverify as av
import math
import re
import subprocess
import json


def _versiontuple(v):
    return tuple(map(int, (v.split("."))))


def evaluate_pyextmod_output(abaverify_obj, jobName, arguments):
    '''
    Helper for evaluating python extension module implementation
    '''
    # Arguments
    subroutine = arguments[0]
    # Run the debug file through the python extension module helper code; outputs a json file
    path_to_debugpy = av.main.options.outputDirectory + '/' + jobName + '-1-debug-0.py'
    subprocess.check_output('cd ../pyextmod && python verify_debug.py ' + subroutine + ' ' + path_to_debugpy, shell=True)
    # Load the json file with the state variables computed by abaqus and the Python Extension Module
    with open(os.path.join(av.main.options.outputDirectory, jobName+'_pyextmod_results.json'), 'r') as f:
        results_dict = json.load(f)
    # Load the file that specifies which state variables to compare (and tolerances for comparison)
    results_expected = __import__('verify_debug_' + jobName + '_expected').parameters
    for sv in results_expected.keys():
        abaverify_obj.assertAlmostEqual(results_dict[sv][0], results_dict[sv][1], delta=results_expected[sv]) # First value is abaqus, 2nd value is pyextmod


def plotFailureEnvelope(baseName, abscissaIdentifier, ordinateIdentifier, abcissaStrengths, ordinateStrengths):
    """
    Create a plot of the failure envelope. Does nothing if matplotlib import fails.
    """

    # Try to import matplotlib
    try:
        import matplotlib as mpl
        if os.environ.get('DISPLAY', '') == '':
            mpl.use('Agg')
        import matplotlib.pyplot as plt

        # Read the failure envelope data
        with open(os.path.join(av.main.options.outputDirectory, baseName + '_failure_envelope.txt'), 'r') as fe:
            data = dict()
            dataHeaders = list()
            for line in fe:
                lineSplit = line.split(', ')

                # Handle the header row separately
                if len(data) == 0:
                    for i in range(0, len(lineSplit)):
                        data[lineSplit[i]] = list()
                        dataHeaders.append(lineSplit[i])
                else:
                    for i in range(0, len(lineSplit)):
                        data[dataHeaders[i]].append(float(lineSplit[i]))

        # Plot the failure envelope
        fig, ax = plt.subplots()

        # Reference data
        dataRef = dict()
        dataRef[abscissaIdentifier] = abcissaStrengths + [0]*len(ordinateStrengths)
        dataRef[ordinateIdentifier] = [0]*len(abcissaStrengths) + ordinateStrengths
        plt.plot(dataRef[abscissaIdentifier], dataRef[ordinateIdentifier], 'x', markeredgecolor='black')

        # Data from CompDam
        plt.plot(data[abscissaIdentifier], data[ordinateIdentifier], 'o', markerfacecolor='none', markeredgecolor='#ED7D31')
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        plt.xlabel(r'$\sigma_{' + abscissaIdentifier.split('S')[1] + '}$ [MPa]')
        plt.ylabel(r'$\sigma_{' + ordinateIdentifier.split('S')[1] + '}$ [MPa]')
        fig.savefig(os.path.join(av.main.options.outputDirectory, baseName + '.png'), dpi=300)

    # If import fails, the above code is skipped
    except ImportError:
        print("INFO: matplotlib package not found. Install matplotlib to generate plots of the failure envelope automatically.")


def plotStressLife(baseName, stressRatios, R_ratio=0.1):
    """
    Create a stress-life plot from a series of fatigue analyses.
    """

    import math
    # Try to import matplotlib
    try:
        import matplotlib as mpl
    # If import fails, the above code is skipped
    except ImportError:
        print("INFO: matplotlib package not found. Install matplotlib to generate stress-life plots automatically.")
        raise
    # if os.environ.get('DISPLAY', '') == '':
    #    mpl.use('Agg')
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    fatigue_life = list()

    for sr in stressRatios:
        # Read the fatigue cycle data from the inc2cycles log files. Assume the last line represents failure.
        file_name = baseName + '_stress_ratio_' + str(sr).replace('.', 'p') + '_inc2cycles.log'
        try:
            with open(os.path.join(av.main.options.outputDirectory, file_name), 'r') as f:
                lines = f.read().splitlines()
                fatigue_life.append(float(lines[-1].split()[-1]))
        except FileNotFoundError:
            print('Unable to find: ', file_name)
    if len(fatigue_life) < len(stressRatios):
        print('Missing results for some of the parametric tests, skipping stress life plot generation')
        return

    # Create the stress-life plot
    fig, ax = plt.subplots()

    abscissaIdentifier = 'Life'
    ordinateIdentifier = 'Fatigue Strength'

    # Analysis data
    data = dict()
    data[abscissaIdentifier] = fatigue_life
    data[ordinateIdentifier] = stressRatios

    plt.plot(data[abscissaIdentifier], data[ordinateIdentifier], 'x', markeredgecolor='black')

    with open(os.path.join(av.main.options.outputDirectory, baseName + '_R-ratio_' + str(R_ratio).replace('.', '') + '.txt'), 'w') as f:
        f.write("{},{}\n".format(ordinateIdentifier, abscissaIdentifier))
        for i in range(len(data[ordinateIdentifier])):
            f.write("{},{}\n".format(data[ordinateIdentifier][i], data[abscissaIdentifier][i]))

    # Formatting
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlim(xmin=1e0, xmax=1e8)
    ax.set_ylim(ymin=0.5, ymax=1.0)
    formatter = FuncFormatter(lambda y, _: '{:.1g}'.format(y))
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.set_minor_formatter(formatter)

    plt.xlabel(abscissaIdentifier + ', cycles')
    plt.ylabel(ordinateIdentifier + ', ' + r'$\sigma^{max}/\sigma_{c}$')
    # plt.grid(b=True, which='major', color='xkcd:silver', linestyle='--')

    fig.savefig(os.path.join(av.main.options.outputDirectory, baseName + '_R-ratio_' + str(R_ratio).replace('.', '') + '.png'), dpi=300)


class ParametricMixedModeMatrix(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    Parametric mixed mode tests.
    """

    # Refers to the template input file name
    baseName = "test_C3D8R_mixedModeMatrix"

    # Range of parameters to test; all combinations are tested
    # alpha is the angle of the crack normal
    # beta defines the direction of tensile loading in Step-1 and compressive loading in Step-2
    parameters = {'alpha': range(0,50,10), 'beta': range(0,210,30), 'friction': [0.00, 0.15, 0.30, 0.45, 0.60]}

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D.inp')


class ParametricElementSizeQuad(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    vucharlength() tests for quad and hex elements not aligned with fiber material direction.
    """

    # Refers to the template input file name
    baseName = "test_C3D8R_elementSize"

    # The angle of misalignment (psi) and the aspect ratio, i.e., L1/L2, (alpha) of the element edges are here varied.
    # A misalignment angle of zero will result in an Abaqus pre error due to an *NMAP rotation command being used in the input deck
    # parameters = {'misalignment_angle': [-45, -30, -15, 1, 11.25, 22.5, 45], 'alpha': [1.0, 1.5]}
    parameters = {'misalignment_angle': [-45,], 'alpha': [1.0, 1.5]}

    # Closed-form equation for the matrix characteristic element length, valid for misalignment angles between -45 and +45 degrees and gamma = 90deg
    L2 = 0.2  # matrix-direction element edge length
    Lc_eq = lambda L, alpha, psi: L * (alpha * math.sin(abs(math.radians(psi))) + math.cos(math.radians(psi)))
    Lc1 = []
    Lc2 = []
    for alpha in parameters['alpha']:
        for psi in parameters['misalignment_angle']:
            Lc1.append(Lc_eq(L2*alpha, 1.0/alpha, psi))
            Lc2.append(Lc_eq(L2, alpha, psi))

    # Element sizes are dependent on the misalignment and skew angles
    expectedpy_parameters = {'Lc1': Lc1, 'Lc2': Lc2}

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D.inp')


class ParametricElementSizeTri(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    vucharlength() tests for tri and wedge elements not aligned with fiber material direction.
    """

    # Refers to the template input file name
    baseName = "test_C3D6_elementSize"

    # The angle of misalignment (psi) and the aspect ratio, i.e., L1/L2, (alpha) of the element edges are here varied.
    # A misalignment angle of zero will result in an Abaqus pre error due to an *NMAP rotation command being used in the input deck
    parameters = {'misalignment_angle': [-45, -30, -15, 1, 5, 10, 15], 'alpha': [1.0, 1.5]}
    # The maximum misalignment_angle value must be less than 0.5*atan(1/alpha) to pass the below test

    # Closed-form equation for the matrix characteristic element length, valid for misalignment angles between -45 and +45 degrees and gamma = 90deg
    L2 = 0.2  # matrix-direction element edge length
    Lc_eq = lambda L, alpha, psi: L * (alpha * math.sin(abs(math.radians(psi))) + math.cos(math.radians(psi)))
    Lc1 = []
    Lc2 = []
    for psi in parameters['misalignment_angle']:
        for alpha in parameters['alpha']:
            Lc1.append(Lc_eq(L2*alpha, 1.0/alpha, psi))
            Lc2.append(Lc_eq(L2, alpha, psi))

    # Element sizes are dependent on the misalignment and skew angles
    expectedpy_parameters = {'Lc1': Lc1, 'Lc2': Lc2}

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D.inp')


class ParametricStressLife(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    Generate data for a stress life plot with a series of fatigue analyses.
    """

    # Refers to the template input file name
    baseName = "test_COH3D8_fatigue_normal"

    parameters = {'stress_ratio': [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]}

    expectedpy_parameters = {'stress_ratio': parameters['stress_ratio']}

    fatigue_R_ratio = 0.5

    @classmethod
    def tearDownClass(cls):
        plotStressLife(baseName=cls.baseName, stressRatios=cls.parameters['stress_ratio'], R_ratio=cls.fatigue_R_ratio)

    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-COH-fatigue.inp')


class ParametricFailureEnvelope_sig12sig22(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    Generate failure envelope in the sigma12 - sigma22 space with a C3D8R element
    """

    # Refers to the template input file name
    baseName = "test_C3D8R_failureEnvelope_sig12sig22"

    # Range of parameters to test; all combinations are tested
    abcissaStrengths = [-199.8, 62.3]
    ordinateStrengths = [92.3]
    parameters = {'loadRatio':  [x/100. for x in range(0,101,5)],
                  'matrixStrength': abcissaStrengths}

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D.inp')

    @classmethod
    def tearDownClass(cls):
        plotFailureEnvelope(baseName=cls.baseName, abscissaIdentifier='S22', ordinateIdentifier='S12', abcissaStrengths=cls.abcissaStrengths, ordinateStrengths=cls.ordinateStrengths)


class ParametricFailureEnvelope_sig11sig22(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    Generate failure envelope in the sigma11 - sigma22 space with C3D8R element
    """

    # Refers to the template input file name
    baseName = "test_C3D8R_failureEnvelope_sig11sig22"

    # Range of parameters to test; all combinations are tested
    ordinateStrengths = [-1200.1, 2326.2]
    abcissaStrengths = [-199.8, 62.3]
    parameters = {'loadRatio':  [x/100. for x in range(0,101,50)],
                  'fiberStrength': ordinateStrengths,
                  'matrixStrength': abcissaStrengths}

    expectedpy_parameters = 'parameters'  # Use same parameters for expected values

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D.inp')

    @classmethod
    def tearDownClass(cls):
        plotFailureEnvelope(baseName=cls.baseName, abscissaIdentifier='S22', ordinateIdentifier='S11', abcissaStrengths=cls.abcissaStrengths, ordinateStrengths=cls.ordinateStrengths)


class ParametricKinkBandWidth_twoElement(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    Tests for fiber compression damage mode to ensure mesh objectivity
    Should yield the same response as ParametricKinkBandWidth_singleElement
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D-25sv.inp')

    # Refers to the template input file name
    baseName = "test_C3D8R_twoElement_fiberCompression_FKT"

    # Use python script instead of input file
    pythonScriptForModel = True

    # Range of parameters to test; all combinations are tested
    parameters = {'elasticElToTotal': [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]}

    # Crush stress is different for each kinkband size, so the expected values are specified here
    expectedpy_parameters = {'crushStress': [-7.9, -8.8, -9.6, -10.3, -11, -11.5]}


class ParametricKinkBandWidth_singleElement(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    Tests to show the effect of kinkband width relative to element size
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D-25sv.inp')

    # Refers to the template input file name
    baseName = "test_C3D8R_fiberCompression_FKT_12"

    # Range of parameters to test; all combinations are tested
    parameters = {'wkbToTotal': [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]}

    # Crush stress is different for each kinkband size, so the expected values are specified here
    expectedpy_parameters = {'crushStress': [-7.9, -8.8, -9.6, -10.3, -11, -11.5, -11.9]}

    # Ignore the displacement at peak load
    expectedpy_ignore = ('x_at_peak_in_xy')


class VerifyDebugPy(av.TestCase):
    """
    Tests to verify that the logic that writes the debug.py file is working
    properly.

    Single element tests are run with a special flag set
    (debug_kill_at_total_time) to trigger an error and write the debug.py
    file. The debug.py file is loaded and executed by the Python Extension
    Module. The tests compare the state variables calculated using Abaqus to
    those calculated using the Python Extension Module. Only the state
    variables listed in the corresponding verify_debug_[jobname]_expected.py
    are checked; the state variables must be equal to pass the test.

    Several single element tests are used to exercise all the state variables.

    The Python Extension Module helper routine 'verify_debug.py' produces a
    json file with the state variables calculated from Abaqus and the Python
    Extension Module. The first entry is the value from Abaqus and the second
    entry is the value from the Python Extension Module.
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles(['IM7-8552-C3D.inp', 'IM7-8552-C3D-25sv.inp'])

        # Check for bash
        try:
            with open(os.path.join(os.getcwd(), os.pardir, 'etc', 'config.json'), 'r') as f:
                config = json.load(f)
                if not config["bash"]:
                    raise av.unittest.SkipTest("CompDam configuration has bash disabled; skipping")
        except IOError:
            print("WARNING: CompDam configuration was not set during installation. Run `python setup.py` from the CompDam root folder.")
            raise av.unittest.SkipTest("WARNING: Bash not found, skipping VerifyDebugPy")

        # Check for abaverify >= 0.5.0
        installed_version = av.__version__
        minimum_required_version = "0.5.0"
        if _versiontuple(installed_version) < _versiontuple(minimum_required_version):
            raise av.unittest.SkipTest("Abaverify 0.5.0 or newer required for these tests; skipping")

        # Run compile script
        sys.stdout.write('Compiling CompDam into a Python Extension Module ... ')
        subprocess.check_output('cd ../pyextmod && make clean && make', stderr=subprocess.STDOUT, shell=True)
        sys.stdout.write(' DONE\n')


    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_C3D8R_matrixTension_debug(self):
        self.runTest("test_C3D8R_matrixTension_debug", func=evaluate_pyextmod_output, arguments=['dgdevolve'],
                     inpName="test_C3D8R_matrixTension",
                     substitutions=[('debug_kill_at_total_time, -1.d0', 'debug_kill_at_total_time, 0.09d0'),])

    def test_C3D8R_fiberCompression_FKT_12_debug(self):
        self.runTest("test_C3D8R_fiberCompression_FKT_12_debug", func=evaluate_pyextmod_output, arguments=['dgdkinkband'],
                     inpName="test_C3D8R_fiberCompression_FKT_12",
                     substitutions=[('debug_kill_at_total_time, -1.d0', 'debug_kill_at_total_time, 0.08d0'),])


class SingleElementSchaeferTests(av.TestCase):
    """
    Single element models to test the CompDam_DGD code base for Schaefer theory responses
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles('IM7-8552-C3D-schaefer.inp')

    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_C3D8R_schaefer_oat90(self):
        """ Simple tension applied in the matrix direction, solid element. Tests Schaefer theory """
        self.runTest("test_C3D8R_schaefer_oat90")

    def test_C3D8R_schaefer_oat30(self):
        """ Off Axis Tension (30 deg) solid element. Tests Schaefer theory """
        self.runTest("test_C3D8R_schaefer_oat30")

    def test_C3D8R_schaefer_oat60(self):
        """ Off Axis Tension (60 deg) solid element. Tests Schaefer theory """
        self.runTest("test_C3D8R_schaefer_oat60")

    def test_C3D8R_schaefer_oat75(self):
        """ Off Axis Tension (75 deg) solid element.  Tests Schaefer theory """
        self.runTest("test_C3D8R_schaefer_oat75")


class SingleElementCohesiveTests(av.TestCase):
    """
    Single element models to test the cohesive element material model features
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles(['test_COH_outputs.inp', 'IM7-8552-COH.inp'])

    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_COH2D4_normal(self):
        """ Single 2-D cohesive element test for normal loading """
        self.runTest("test_COH2D4_normal")

    def test_COH2D4_shear(self):
        """ Single 2-D cohesive element test for shear loading """
        self.runTest("test_COH2D4_shear")

    def test_COH2D4_shear_compression(self):
        """ Single 2-D cohesive element test for shear loading with normal compression """
        self.runTest("test_COH2D4_shear_compression")

    def test_COH2D4_shear_friction(self):
        """ Single 2-D cohesive element test for shear loading with friction """
        self.runTest("test_COH2D4_shear_friction")

    def test_COH3D8_normal(self):
        """ Single cohesive element test for mode I response """
        self.runTest("test_COH3D8_normal")

    def test_COH3D8_normal_deletion(self):
        """ Single cohesive element test for mode I response with element deletion """
        self.runTest("test_COH3D8_normal_deletion")

    def test_COH3D8_normal_constitutive_thickness(self):
        """ Single cohesive element test for mode I response with nondefault constitutive thickness """
        self.runTest("test_COH3D8_normal_constitutive_thickness",
                     inpName="test_COH3D8_normal_effmod",
                     substitutions=[('constit_thickness = 1.0', 'constit_thickness = calc_ct'),
                                    ('density = mat_density', 'density = calc_rho'),
                                    ('mass_scaling_factor = 1000', 'mass_scaling_factor = 1')],
                     expected_substitutions=[('(0.03, 1.05e-4)', '(0.001, 2e-8)')  # DT assertion consistent with target_dt and mass_scaling_factor = 1
                                            ])

    def test_COH3D8_normal_effmod(self):
        """ Single cohesive element test for mode I response using EFFMOD """
        self.runTest("test_COH3D8_normal_effmod")

    def test_COH3D8_normal_effmod_ct_0p2(self):
        """ Single cohesive element test for mode I response with constitutive thickness = 0.2 """
        self.runTest("test_COH3D8_normal_effmod_ct_0p2",
                     inpName="test_COH3D8_normal_effmod",
                     substitutions=[('constit_thickness = 1.0', 'constit_thickness = 0.2'),],
                     datacheck=True)

    def test_COH3D8_normal_effmod_ct(self):
        """ Single cohesive element test for mode I response with target dt = 2e-8 """
        self.runTest("test_COH3D8_normal_effmod_ct",
                     inpName="test_COH3D8_normal_effmod",
                     substitutions=[('constit_thickness = 1.0', 'constit_thickness = calc_ct'),
                                    ('density = mat_density', 'density = calc_rho'),],
                     datacheck=True)

    def test_COH3D8_normal_effmod_ct_halfYT(self):
        """ Single cohesive element test for mode I response with target dt = 2e-8, using 50% YT """
        self.runTest("test_COH3D8_normal_effmod_ct_halfYT",
                     inpName="test_COH3D8_normal_effmod",
                     substitutions=[('YT = 62.3', 'YT = 31.15'),
                                    ('constit_thickness = 1.0', 'constit_thickness = calc_ct'),
                                    ('density = mat_density', 'density = calc_rho'),],
                     datacheck=True)

    def test_COH3D8_shear13(self):
        """ Single cohesive element test for 1-3 shear loading """
        self.runTest("test_COH3D8_shear13")

    def test_COH3D8_shear13_effmod(self):
        """ Single cohesive element test for 1-3 shear loading using EFFMOD """
        self.runTest("test_COH3D8_shear13_effmod",
                     inpName="test_COH3D8_shear13",
                     substitutions=[('material_options = ""', 'material_options = ", effmod"'),],
                     expected_substitutions=[('\"tolerance\": 0.2', '"tolerance": 1.0'),
                                             ('disp_at_zero_tol_pct = 0.001', 'disp_at_zero_tol_pct = 0.01')])

    def test_COH3D8_shear13_constitutive_thickness(self):
        """ Single cohesive element test for 1-3 shear loading with nondefault constitutive thickness """
        self.runTest("test_COH3D8_shear13_constitutive_thickness",
                     inpName="test_COH3D8_shear13",
                     substitutions=[('material_options = ""', 'material_options = ", effmod"'),
                                    ('constit_thickness = 1.0', 'constit_thickness = 74.59'),
                                    ('density = 1.57e-09', 'density = 4.21e-13')],
                     expected_substitutions=[('\"tolerance\": 0.2', '"tolerance": 1.0'),
                                             ('disp_at_zero_tol_pct = 0.001', 'disp_at_zero_tol_pct = 0.01')])

    def test_COH3D8_shear13_compression(self):
        """ Single cohesive element test for 1-3 shear loading with normal compression """
        self.runTest("test_COH3D8_shear13_compression")

    def test_COH3D8_shear13_compression_effmod(self):
        """ Single cohesive element test for 1-3 shear loading with normal compression using EFFMOD """
        self.runTest("test_COH3D8_shear13_compression_effmod",
                     inpName="test_COH3D8_shear13_compression",
                     substitutions=[('material_options = ""', 'material_options = ", effmod"'),],
                     expected_substitutions=[('\"tolerance\": 0.2', '"tolerance": 1.0'),
                                             ('disp_at_zero_tol_pct = 0.001', 'disp_at_zero_tol_pct = 0.01')])

    def test_COH3D8_shear13_friction(self):
        """ Single cohesive element test for 1-3 shear loading with friction """
        self.runTest("test_COH3D8_shear13_friction")

    def test_COH3D8_shear23(self):
        """ Single cohesive element test for 2-3 shear loading """
        self.runTest("test_COH3D8_shear23")

    def test_COH3D8_shear23_compression(self):
        """ Single cohesive element test for 2-3 shear loading with normal compression """
        self.runTest("test_COH3D8_shear23_compression")

    def test_COH3D8_shear23_friction(self):
        """ Single cohesive element test for 2-3 shear loading with friction """
        self.runTest("test_COH3D8_shear23_friction")
        
    def test_COH3D8_thick_normal(self):
        """ Single cohesive element test for finite-thickness mode I response """
        self.runTest("test_COH3D8_thick_normal")
        
    def test_COH3D8_thick_shear13(self):
        """ Single cohesive element test for finite-thickness 1-3 shear loading """
        self.runTest("test_COH3D8_thick_shear13")
        
    def test_COH3D8_thick_shear23(self):
        """ Single cohesive element test for finite-thickness 2-3 shear loading """
        self.runTest("test_COH3D8_thick_shear23")


class SingleElementTests(av.TestCase):
    """
    Single element models to test the solid element material model features
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles(['IM7-8552-C3D.inp', 'IM7-8552-C3D-25sv.inp', 'IM7-8552-C3D-schaefer.inp',
                      '_default_parameters.inp'])

    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_C3D8R_error(self):
        """ Intentionally cause a DGD convergence error """
        self.runTest("test_C3D8R_error")


    def test_C3D8R_matrixTension(self):
        """ Simple tension in the matrix direction, with damage """
        self.runTest("test_C3D8R_matrixTension")


    def test_C3D8R_matrixTension_deletion(self):
        """ Simple tension in the matrix direction, with damage and element deletion """
        self.runTest("test_C3D8R_matrixTension_deletion")


    def test_C3D8R_matrixTensionFMS(self):
        """ Simple tension in the matrix direction, with damage, large opening, fixed mass scaling """
        self.runTest("test_C3D8R_matrixTensionFMS")


    def test_C3D8R_matrixTensionVMS(self):
        """ Simple tension in the matrix direction, with damage, large opening, variable mass scaling """
        self.runTest("test_C3D8R_matrixTensionVMS",
                     inpName="test_C3D8R_matrixTensionFMS",
                     substitutions=[('Fixed mass scaling', 'Variable mass scaling, freq=1')],
                     expected_substitutions=[('\"tolerance\": 0.5', '\"tolerance\": 0.6')],  # Tolerance for S22 continuous
                     )


    def test_C3D8R_matrixTension90(self):
        """ Simple tension in the 3-direction, with damage """
        self.runTest("test_C3D8R_matrixTension90")


    @av.unittest.expectedFailure
    def test_C3D8R_matrixTension90_missingIC(self):
        """ Simple tension in the 3-direction, with damage, missing initial conditions """
        self.runTest("test_C3D8R_matrixTension90_missingIC",
                     inpName="test_C3D8R_matrixTension90",
                     substitutions=[("\*Initial Conditions, Type=Solution", "** *Initial Conditions, Type=Solution"),
                                    (" ALL_ELEMS,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,", "**  ALL_ELEMS,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,"),
                                    (" 0.d0,       0.d0, 90.d0,     1,  0.d0,  0.d0,  0.d0,  0.d0,", "**  0.d0,       0.d0, 90.d0,     1,  0.d0,  0.d0,  0.d0,  0.d0,"),
                                    (" 0.d0,       0.d0,  0.d0,  0.d0", "**  0.d0,       0.d0,  0.d0,  0.d0")]
                    )


    def test_C3D8R_simpleShear12(self):
        """ Simple shear in the 1-2 plane, with damage """
        self.runTest("test_C3D8R_simpleShear12")


    def test_C3D8R_simpleShear12friction(self):
        """ Compression followed by simple shear in the 1-2 plane """
        self.runTest("test_C3D8R_simpleShear12friction")


    def test_C3D8R_nonlinearShear12(self):
        """ Nonlinear shear model, loading and unloading """
        self.runTest("test_C3D8R_nonlinearShear12")


    def test_C3D8R_nonlinearShear12_monotonic(self):
        """ Nonlinear shear model, monotonic loading """
        self.runTest("test_C3D8R_nonlinearShear12_monotonic")


    def test_C3D8R_nonlinearShear12_withFKT(self):
        """ Nonlinear shear model, loading and unloading, with FKT """
        self.runTest("test_C3D8R_nonlinearShear12_withFKT")


    def test_C3D8R_nonlinearShear12_loadReversal(self):
        """ Nonlinear shear model, loading and unloading, no damage, including full load reversal """
        self.runTest("test_C3D8R_nonlinearShear12_loadReversal")


    def test_C3D8R_fiberTension(self):
        """ Simple tension in fiber direction, with damage """
        self.runTest("test_C3D8R_fiberTension")


    def test_C3D8R_fiberTension_deletion(self):
        """ Simple tension in fiber direction, with damage and element deletion """
        self.runTest("test_C3D8R_fiberTension_deletion")


    def test_C3D8R_fiberTension_FN(self):
        """ Simple tension in fiber direction, with damage and fiber nonlinearity """
        self.runTest("test_C3D8R_fiberTension_FN")


    def test_C3D8R_fiberCompression_FKT_12(self):
        """ Fiber compression: Fiber kinking theory based model, 1-2 """
        self.runTest("test_C3D8R_fiberCompression_FKT_12")


    def test_C3D8R_fiberCompression_FKT_13(self):
        """ Fiber compression: Fiber kinking theory based model, 1-3 """
        self.runTest("test_C3D8R_fiberCompression_FKT_13")


    def test_C3D8R_fiberCompression_FKT_3D(self):
        """ Fiber compression: Fiber kinking theory based model, 3-D """
        self.runTest("test_C3D8R_fiberCompression_FKT_3D")


    def test_C3D8R_fiberCompression_FKT_12_FF(self):
        """ Fiber compression: Fiber kinking theory based model, fiber failure """
        self.runTest("test_C3D8R_fiberCompression_FKT_12_FF",
                     inpName="test_C3D8R_fiberCompression_FKT_12",
                     substitutions=[('fkt_fiber_failure_angle, -1.d0', 'fkt_fiber_failure_angle, 10.d0')])


    def test_C3D8R_fiberCompression_FKT_13_FF(self):
        """ Fiber compression: Fiber kinking theory based model, fiber failure """
        self.runTest("test_C3D8R_fiberCompression_FKT_13_FF",
                     inpName="test_C3D8R_fiberCompression_FKT_13",
                     substitutions=[('fkt_fiber_failure_angle, -1.d0', 'fkt_fiber_failure_angle, 10.d0')])

    def test_C3D8R_fiberCompression_FKT_12_FF_negphi0(self):
        """ Fiber compression: Fiber kinking theory based model, fiber failure """
        self.runTest("test_C3D8R_fiberCompression_FKT_12_FF_negphi0")


    def test_C3D8R_fiberCompression_FKT_12_FN(self):
        """ Fiber compression: Fiber kinking theory based model, fiber nonlinearity """
        self.runTest("test_C3D8R_fiberCompression_FKT_12_FN")


    def test_C3D8R_fiberCompression_FKT_13_FN(self):
        """ Fiber compression: Fiber kinking theory based model, fiber nonlinearity """
        self.runTest("test_C3D8R_fiberCompression_FKT_13_FN")


    def test_C3D8R_fiberCompression_FKT_3D_pert(self):
        """ Fiber compression: Fiber kinking theory based model, 3-D, perturbation """
        self.runTest("test_C3D8R_fiberCompression_FKT_3D_pert")


    def test_C3D8R_fiberCompression_FKT_3D_spring_oop(self):
        """ Fiber compression: Fiber kinking theory based model, 3-D, out-of-plane stiffness """
        self.runTest("test_C3D8R_fiberCompression_FKT_3D_spring_oop")


    def test_C3D8R_fiberCompression_FKT_3D_spring_ip(self):
        """ Fiber compression: Fiber kinking theory based model, 3-D, in-plane stiffness """
        self.runTest("test_C3D8R_fiberCompression_FKT_3D_spring_ip")


    def test_C3D8R_fiberCompression_BL(self):
        """ Fiber compression: Bilinear softening based model """
        self.runTest("test_C3D8R_fiberCompression_BL")


    def test_C3D8R_fiberCompression_BL_deletion(self):
        """ Fiber compression: Bilinear softening based model, with element deletion """
        self.runTest("test_C3D8R_fiberCompression_BL_deletion")


    def test_C3D8R_fiberCompression_BL_FN(self):
        """ Fiber compression: Bilinear softening based model, fiber nonlinearity """
        self.runTest("test_C3D8R_fiberCompression_BL_FN")


    def test_C3D8R_fiberLoadReversal(self):
        """ Fiber damage model, Maimi: load reversal """
        self.runTest("test_C3D8R_fiberLoadReversal")


    def test_C3D8R_fiberLoadReversal_FN(self):
        """ Fiber damage model, Maimi: load reversal, fiber nonlinearity """
        self.runTest("test_C3D8R_fiberLoadReversal_FN")


    def test_C3D8R_nonlinearShear12(self):
        """ Nonlinear shear model, loading and unloading in 1-2 plane """
        self.runTest("test_C3D8R_nonlinearShear12")


    def test_C3D8R_nonlinearShear12_loadReversal(self):
        """ Nonlinear shear model, loading and unloading in 1-2 plane, including full load reversal """
        self.runTest("test_C3D8R_nonlinearShear12_loadReversal")


    def test_C3D8R_nonlinearShear13(self):
        """ Nonlinear shear model, loading and unloading in 1-3 plane"""
        self.runTest("test_C3D8R_nonlinearShear13")


    def test_C3D8R_nonlinearShear13_loadReversal(self):
        """ Nonlinear shear model, loading and unloading in 1-3 plane, including full load reversal"""
        self.runTest("test_C3D8R_nonlinearShear13_loadReversal")


    def test_C3D8R_schapery12(self):
        """ Schapery micro-damage model, loading and unloading in 1-2 plane"""
        self.runTest("test_C3D8R_schapery12")


    def test_C3D8R_matrixCompression(self):
        """ Simple compression in the matrix direction """
        self.runTest("test_C3D8R_matrixCompression")
        
        
    def test_C3D8R_matrixCompression_friction(self):
        """ Simple compression in the matrix direction with friction"""
        self.runTest("test_C3D8R_matrixCompression_friction")


    def test_C3D8R_elastic_matrixTension(self):
        """ Simple tension in the matrix direction, no damage """
        self.runTest("test_C3D8R_elastic_matrixTension")


    def test_C3D8R_elastic_fiberTension(self):
        """ Simple tension in the fiber direction, no damage """
        self.runTest("test_C3D8R_elastic_fiberTension")


    def test_C3D8R_elastic_simpleShear12(self):
        """ Simple shear in the 1-2 plane, no damage """
        self.runTest("test_C3D8R_elastic_simpleShear12")


    def test_CPS4R_elementSize(self):
        """ Characteristic element size test, plane stress element """
        self.runTest("test_CPS4R_elementSize")


    def test_C3D6_matrixTension(self):
        """ Simple tension in the matrix direction, two wedge elements """
        self.runTest("test_C3D6_matrixTension")


    def test_C3D6_simpleShear12(self):
        """ Simple shear in the 1-2 plane, two wedge elements """
        self.runTest("test_C3D6_simpleShear12")


    def test_S4R_elementSize(self):
        """ Characteristic element size test, conventional shell element """
        self.runTest("test_S4R_elementSize")
        
        
    def test_C3D8R_residualStress(self):
        """ Residual thermal stress in a solid element """
        self.runTest("test_C3D8R_residualStress")


class SingleElementFatigueTests(av.TestCase):
    """
    Single element models to test the matrix crack fatigue model
    """
    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_COH3D8_fatigue_normal(self):
        """ Single cohesive element fatigue test for mode I loading """
        self.runTest("test_COH3D8_fatigue_normal")

    def test_COH3D8_fatigue_shear13(self):
        """ Single cohesive element fatigue test for 1-3 shear loading """
        self.runTest("test_COH3D8_fatigue_shear13")
        
    def test_C3D8R_fatigue_matrixTension(self):
        """ Single solid element fatigue test for tensile matrix loading """
        self.runTest("test_C3D8R_fatigue_matrixTension")

    def test_C3D8R_fatigue_simpleShear12(self):
        """ Single solid element fatigue test for simple shear loading in the 1--2 plane """
        self.runTest("test_C3D8R_fatigue_simpleShear12")


if __name__ == "__main__":
    av.runTests(relPathToUserSub='../for/CompDam_DGD', double=True)
