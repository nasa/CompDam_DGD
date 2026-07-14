#
# Code to run double cantilever beam (DCB) example problems
#

import os
import abaverify as av
import textwrap

include_files = ["_DCB_3D_parameters.inp", "_DCB_3D_geometry+sec.inp", "_DCB_3D_exp_step.inp",
                 "_outputs.inp", "_outputs_compdam.inp", "_interface_C3D8R.inp", "_interface_COH3D8.inp",
                 "_interface_COH3D8_insert.inp", "_interface_tie.inp", "_IM7_8552.inp",]

def parameter_keyword_to_dict(inpfile):
    """
    Reads parameter keywords in the provided inpfile and returns a dictionary
    """
    with open(inpfile) as f:
        lines = f.readlines()
    parameters = {}
    is_parameter_dataline = False
    for line in lines:
        if line.startswith('*Parameter'):
            is_parameter_dataline = True
            continue
        elif line.startswith('**'):
            continue
        elif line.startswith('*'):
            is_parameter_dataline = False
        if is_parameter_dataline:
            ls = line.split('=')
            key = ls[0].strip()
            value = ls[1].strip()
            if '#' in value:
                value = value.split('#')[0].strip()
            try:
                value = float(value)
            except:
                pass
            parameters[key] = value
    return parameters

def make_dynamic():
    av.copyFiles(
        ["_DCB_3D_parameters.inp", "_DCB_3D_exp_step.inp", "_IM7_8552.inp"], 
        substitutions=[
            [
                ('TL = 120.0', 'TL = 25.0'),
                ('crackLength = 50.8', 'crackLength = 10.0'),
                ('FM_after = 65.0', 'FM_after = 12.0'),
                ('arm_num_plies = 19', 'arm_num_plies = 1'),
                ('disp = 4.5', 'disp = 6.0'),
                ('step_time = 1', 'step_time = 0.002'),
                ('ramp_time_frac = 0.1', 'ramp_time_frac = 0.3'),
                ('number_frames = 50', 'number_frames = 100')
            ],[
                ('T_Node, 3, 3, 0.5', 'T_Node, 3, 3, 1.0'), 
                ('B_Node, 3, 3, -0.5', 'B_Node, 3, 3, 0.'),
            ],[
                ('YT    = 62.3', 'YT    = 30.0'),
                ('SL    = 92.30', 'SL    = 30.0'),
                ('Pen2    = 2.e5', 'Pen2    = 1.e5'),
            ]
        ]
    )

def make_finite_thickness(interface_el_thk, param, include_element_deletion=False):
    nom_arm_thk = param['ply_thk']*param['arm_num_plies']
    # scaled material properties
    sf = nom_arm_thk/(nom_arm_thk-interface_el_thk/2)
    E11 = param['E11']*sf
    E22 = param['E22']*sf
    E33 = E22
    G12 = param['G12']*sf
    G13 = G12
    G23 = param['E22']/2./(1 + param['nu23'])*sf
    pen_N = param['E22']/interface_el_thk
    pen_shear_L = param['G12']/interface_el_thk
    pen_shear_T = G23/sf/interface_el_thk
    interface_coh3d8_subs = [
                ('\*User Material, constants=40, effmod', '*User Material, constants=40'),  # effmod does nothing for finite thickness cohesive elements
                (' <GIIc>,      <BK>,    <YC>,    <alpha0>,    <Pen2>,        ,         ,          ,',
                    f' <GIIc>,      <BK>,    <YC>,    <alpha0>, {pen_N}, {pen_shear_L}, {pen_shear_T},          ,'), 
            ]
    if include_element_deletion:
        element_deletion_keywords = textwrap.dedent(
            """\
            ** CompDam parameters
            **
            *PARAMETER TABLE TYPE, NAME=KEY_INT, parameters=2
            STRING
            INTEGER
            *TABLE COLLECTION, NAME=COMPDAM_PARAM
            ** CompDam parameters with integer values
            ** 0 = False, 1 = True
            *PARAMETER TABLE, type=KEY_INT, label=CDP_INT
            set_status_0_on_d2, 1
            """
        )
        interface_coh3d8_subs.append(('\*Depvar', '*Depvar, delete=11'),)
        interface_coh3d8_subs.append(('\*Material, name=cohesive\n', element_deletion_keywords+'*Material, name=cohesive\n'),)
    av.copyFiles(
        ["_DCB_3D_geometry+sec.inp", "_interface_COH3D8.inp"], 
        substitutions=[
            [
                # Scaled material properties for arms
                ('<E11>, <E22>, <E33>, <nu12>, <nu13>, <nu23>, <G12>, <G13>',
                    f'{E11}, {E22}, {E33}, <nu12>, <nu13>, <nu23>, {G12}, {G13}'),
                ('<G23>,', f'{G23},'),
            ],interface_coh3d8_subs,
        ],
        apply_to_outputDirectory=True
    )

def make_multimaterial_zones(num_zones=4, crackLength=None, mesh_size_coarse=None, mesh_size_fine=None,
                             FM_before=None, FM_after=None, NEW=None, NET=None, TL=None):
    """
    Create elsets and material cards for cohesive multi-material zones
    """
    #
    # ELSETS
    #
    param = parameter_keyword_to_dict("_DCB_3D_parameters.inp")
    if crackLength is None:
        crackLength = param['crackLength']
    if mesh_size_coarse is None:
        mesh_size_coarse = param['mesh_size_coarse']
    if mesh_size_fine is None:
        mesh_size_fine = param['mesh_size_fine']
    if FM_before is None:
        FM_before = param['FM_before']
    if FM_after is None:
        FM_after = param['FM_after']
    if NEW is None:
        NEW = int(param['NEW'])
    if NET is None:
        NET = int(param['NET'])
    if TL is None:
        TL = param['TL']

    NEA = max(5, int((crackLength - FM_before)/mesh_size_coarse))       # Between loading application and start of fine mesh
    NEB = max(5, int(FM_before/mesh_size_fine))  # In fine mesh before init crack tip
    NEC = max(5, int(FM_after/mesh_size_fine))  # In fine mesh after init crack tip
    NED = max(5, int((TL-(crackLength + FM_after))/mesh_size_coarse))
    NEL = NEA + NEB + NEC + NED
    coh_f = NEL*NEW*NET+1+NEA+NEB
    zone_num_el_x = int(NEC/num_zones)
    element_numbers = [
        dict(start=coh_f,                  end=coh_f+zone_num_el_x,   incrementx=1, incrementy=NEL),
        dict(start=coh_f+zone_num_el_x+1,  end=coh_f+zone_num_el_x*2, incrementx=1, incrementy=NEL),
        dict(start=coh_f+zone_num_el_x*2+1,end=coh_f+zone_num_el_x*3, incrementx=1, incrementy=NEL),
        dict(start=coh_f+zone_num_el_x*3+1,end=coh_f+NEC-1,  incrementx=1, incrementy=NEL),
    ]

    inp_elsets = ""
    for i in range(len(element_numbers)):
        inp_elsets += f"*Elset, elset=CMM_Z{i+1}, generate\n"
        for k in range(NEW):
            start = element_numbers[i]['start'] + k*element_numbers[i]['incrementy']
            end = element_numbers[i]['end'] + k*element_numbers[i]['incrementy']
            inp_elsets += f" {start}, {end}\n"

    # Remove last newline
    inp_elsets = inp_elsets[:-1]

    #
    # MATERIAL CARDS
    #
    param = parameter_keyword_to_dict("_IM7_8552.inp")
    YT = param['YT']
    SL = param['SL']
    GIc = param['GIc']
    GIIc = param['GIIc']
    default_material_defintion = """
    *Cohesive section, elset=COH, response=TRACTION SEPARATION, material=cohesive, controls=EC-1, orientation=OID_COH
    <constitutive_thk>
    *Material, name=cohesive
    *Density
    <density_interface>,
    *Depvar
    19,
    1,  COH_dmg
    2,  COH_delta_s1
    3,  COH_delta_n
    4,  COH_delta_s2
    5,  COH_B
    9,  COH_FI
    11, COH_STATUS
    15, COH_slide1
    16, COH_slide2
    *User Material, constants=40, effmod
    ** feature flags, ,   thickness,          4, 5, 6, 7, 8
            200000, ,   <constitutive_thk>,  ,  ,  ,  ,  ,
    **  9         10        11        12        13        14        15        16
    **  E1,       E2,       G12,      nu12,     nu23,     YT,       SL        GYT
        ,         ,          ,          ,         ,     <YT>,    <SL>,    <GIc>,
    **  17        18        19        20        21        22        23        24
    **  GSL,      eta_BK,   YC,       alpha0    E3,       G13,      G23,      nu13,
    <GIIc>,      <BK>,    <YC>,    <alpha0>,    <Pen2>,        ,         ,          ,
    **  25        26        27        28        29        30        31        32
    **  alpha11,  alpha22,  alpha_PL, n_PL,     XT,       fXT,      GXT,      fGXT,
            ,         ,          ,     ,       ,          ,         ,          ,
    **  33        34        35        36        37        38        39        40
    **  XC,       fXC,      GXC,      fGXC,     phi_ff,   w_kb,     phi0,     mu
        ,          ,         ,          ,           ,       ,         ,     0.0"""
    mat_cards = ""
    factors = [1.2, 1.3, 1.4]
    for i in range(1, num_zones):
        mat_card = default_material_defintion
        mat_card = mat_card.replace('elset=COH', f'elset=CMM_Z{i+1}')
        mat_card = mat_card.replace('material=cohesive', f'material=CMM_Z{i+1}')
        mat_card = mat_card.replace('name=cohesive', f'name=CMM_Z{i+1}')
        mat_card = mat_card.replace('<YT>', f'{YT + factors[i-1]:.2f}')
        mat_card = mat_card.replace('<SL>', f'{SL + factors[i-1]:.2f}')
        mat_card = mat_card.replace('<GIc>', f'{GIc * factors[i-1]**2:.2f}')
        mat_card = mat_card.replace('<GIIc>', f'{GIIc * factors[i-1]**2:.2f}')
        mat_cards += mat_card

    # Add cohesive section for back region
    mat_cards += """
                *Cohesive section, elset=BACK_COH, response=TRACTION SEPARATION, material=cohesive, controls=EC-1, orientation=OID_COH
                1.0"""

    return inp_elsets, textwrap.dedent(mat_cards)


class DCBQuasiStatic(av.TestCase):
    """
    Quasi-static DCB models
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles(include_files)
        cls.param = parameter_keyword_to_dict("_IM7_8552.inp")
        cls.param.update(parameter_keyword_to_dict("test_DCB_3D_compdam.inp"))
        cls.param.update(parameter_keyword_to_dict("_DCB_3D_parameters.inp"))

    @classmethod
    def tearDown(cls):
        av.copyFiles(include_files)

    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_DCB_3D_compdam(self):
        """ Baseline 3D DCB example """
        self.runTest("test_DCB_3D_compdam")
    
    def test_DCB_3D_compdam_23(self):
        """ 3D DCB example with propagation in the cohesive 2-3 orientation """
        av.copyFiles(
            ["_DCB_3D_parameters.inp",],
            substitutions=[[
                ('coh_rot_z = 0', 'coh_rot_z = 90'),
            ],]
        )
        self.runTest("test_DCB_3D_compdam_23",
                     inpName="test_DCB_3D_compdam",)

    def test_DCB_3D_compdam_C3D8R(self):
        """ 3D DCB example, using solid elements for fracture interface """
        self.runTest("test_DCB_3D_compdam_C3D8R",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[("interface_el_type = 'COH3D8'", "interface_el_type = 'C3D8R'"), 
                                    ('interface_el_thk = 0.0', 'interface_el_thk = 0.05'), 
                                    ("sdv_outputs = 'SDV1, SDV2, SDV3, SDV4, SDV5, SDV9, SDV11'", "sdv_outputs = 'SDV'"),
                                    ("input=_interface_COH3D8.inp", "input=_interface_C3D8R.inp"),

                                ],)

    def test_DCB_3D_compdam_finite_thk(self):
        """ 3D DCB example, using finite thickness cohesive elements """
        interface_el_thk = 0.02
        density_coh = self.param['density']*interface_el_thk
        make_finite_thickness(interface_el_thk, self.param, include_element_deletion=True) # Note: element deletion is required to get correct response with extended damage progression
        self.runTest("test_DCB_3D_compdam_finite_thk",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[('interface_el_thk = 0.0', 'interface_el_thk = {}'.format(interface_el_thk)), 
                                    ('density_interface = density', 'density_interface = {}'.format(density_coh)),
                                    ('Variable mass scaling, dt=<time_inc>, type=below min, freq=1', 'Fixed mass scaling, dt=<time_inc>, type=below min')
                                ],)
        av.copyFiles(["_DCB_3D_geometry+sec.inp", "_interface_COH3D8.inp"])


class DCBQuasiStaticInsert(av.TestCase):
    """
    Quasi-static DCB models that include damage initialization and/or element deletion
    for modeling the pre-implanted insert (crack starter)
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles(include_files)

    @classmethod
    def tearDown(cls):
        av.copyFiles(include_files)

    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_DCB_3D_compdam_insert_dmg1(self):
        """ Include cohesive elements in insert region, initialized to fully damaged """
        av.copyFiles(
            ["_outputs_compdam.inp",],
            substitutions=[[
                ("\*\* \*Element Output, elset=COH_INSERT", "*Element Output, elset=COH_INSERT"),
                ("\*\* <sdv_outputs>", "<sdv_outputs>"),
            ],]
        )
        self.runTest("test_DCB_3D_compdam_insert",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[("input=_interface_COH3D8.inp", "input=_interface_COH3D8.inp\n*Include, input=_interface_COH3D8_insert.inp"),
                                ],)

    def test_DCB_3D_compdam_insert_delete(self):
        """ Include cohesive elements in insert region, initialized to fully damaged, element deletion """
        av.copyFiles(
            ["_interface_COH3D8.inp", "_interface_COH3D8_insert.inp","_outputs_compdam.inp",],
            substitutions=[[
                ("\*Depvar", "*Depvar, delete=11")
            ],[
                ("\*Depvar", "*Depvar, delete=11")
            ],[
                ("\*\* \*Element Output, elset=COH_INSERT", "*Element Output, elset=COH_INSERT"),
                ("\*\* <sdv_outputs>", "<sdv_outputs>"),
            ],]
        )
        element_deletion_keywords = textwrap.dedent(
            """\
            ** CompDam parameters
            **
            *PARAMETER TABLE TYPE, NAME=KEY_INT, parameters=2
            STRING
            INTEGER
            *TABLE COLLECTION, NAME=COMPDAM_PARAM
            ** CompDam parameters with integer values
            ** 0 = False, 1 = True
            *PARAMETER TABLE, type=KEY_INT, label=CDP_INT
            set_status_0_on_d2, 1
            """
        )
        self.runTest("test_DCB_3D_compdam_insert_delete",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[("input=_interface_COH3D8.inp", "input=_interface_COH3D8.inp\n*Include, input=_interface_COH3D8_insert.inp"),
                                    ("Include, input=_IM7_8552.inp", "Include, input=_IM7_8552.inp\n"+element_deletion_keywords,),],
                    )


class DCBLarge(av.TestCase):
    """
    DCB models that are large (computationally)
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles(include_files)
        cls.param = parameter_keyword_to_dict("_IM7_8552.inp")
        cls.param.update(parameter_keyword_to_dict("test_DCB_3D_compdam.inp"))
        cls.param.update(parameter_keyword_to_dict("_DCB_3D_parameters.inp"))

    @classmethod
    def tearDown(cls):
        av.copyFiles(include_files)

    # -----------------------------------------------------------------------------------------
    # Test methods
    def _model_setup(self, updated_parameters, width):
        additional_elsets, mat_cards = make_multimaterial_zones(**updated_parameters)
        av.copyFiles(
            ["_outputs_compdam.inp", "_DCB_3D_parameters.inp", "_interface_COH3D8.inp", "_DCB_3D_exp_step.inp"],
            substitutions=[[
                ("\*\* \*Element Output, elset=COH_INSERT", "*Element Output, elset=COH_INSERT"),
                ("\*\* <sdv_outputs>", "<sdv_outputs>"),
            ],[
                ("TL = 120.0", f"TL = {updated_parameters['TL']}"),
                ("FM_after = 65.0", f"FM_after = {updated_parameters['FM_after']}"),
                ("width = 1", f"width = {width}"),
                ("disp = 4.5", "disp = 12"),
                ("NEW = 2", f"NEW = {updated_parameters['NEW']}"),
            ],[
                ("\*Elset, elset=first_coh_element", additional_elsets+"\n*Elset, elset=first_coh_element"),
                ("Cohesive section, elset=COH", "Cohesive section, elset=CMM_Z1")
            ],[
                (" TOP_LAYER, 2, 2, 0.", "** TOP_LAYER, 2, 2, 0."),
                (" BOTTOM_LAYER, 2, 2, 0.", "** BOTTOM_LAYER, 2, 2, 0."),
            ]],
            append=["","",mat_cards,""]
        )

    def test_DCB_3D_compdam_multimaterial_narrow(self):
        """ 3D DCB example, using multiple user materials """
        updated_parameters = {
            'TL': 160.0,
            'FM_after': 101.6,
            'NEW': 2,  # Number of elements across the width
        }
        width = 1.0
        self._model_setup(updated_parameters, width)
        self.runTest(f"test_DCB_3D_compdam_multimaterial_narrow", inpName="test_DCB_3D_compdam_multimaterial",
                     expectedpyName="test_DCB_3D_compdam")

    def test_DCB_3D_compdam_multimaterial_wide(self):
        """ 3D DCB example, using multiple user materials """
        updated_parameters = {
            'TL': 160.0,
            'FM_after': 101.6,
            'NEW': 40,  # Number of elements across the width
        }
        width = 20.0
        self._model_setup(updated_parameters, width)
        options = av.getOptions()
        self.runTest(f"test_DCB_3D_compdam_multimaterial_cpus{options.cpus}", inpName="test_DCB_3D_compdam_multimaterial",
                     expectedpyName="test_DCB_3D_compdam", expected_substitutions=[
                         ("peakLoad = 8.0", f"peakLoad = {8.0*width}"),
                         ("finalLoad = 4.38", f"finalLoad = {4.38*width}")])


class DCBDynamic(av.TestCase):
    """
    Dynamic DCB models
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        av.copyFiles(include_files)
        cls.param = parameter_keyword_to_dict("_IM7_8552.inp")
        cls.param.update(parameter_keyword_to_dict("_DCB_3D_parameters.inp"))
        cls.param.update(parameter_keyword_to_dict("test_DCB_3D_compdam.inp"))

    @classmethod
    def tearDown(cls):
        av.copyFiles(include_files)

    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_DCB_3D_dyn_tie(self):
        """ Dynamic DCB, cohesive interface replaced with a tie constraint"""
        make_dynamic()
        self.runTest("test_DCB_3D_dyn_tie",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[('_interface_COH3D8.inp', '_interface_tie.inp'),
                                    ('\*Include, input=_outputs_compdam.inp', '**'),
                                    ('\*Variable mass scaling.*', '**'),],   # No mass scaling
                     expectedpyName="test_DCB_3D_dynamic"
                    )

    @av.unittest.expectedFailure  # This case has significant mass in the coh element, and a different dynamic response
    def test_DCB_3D_dyn_coh(self):
        """ Dynamic DCB, regular zero-thickness cohesive interface"""
        make_dynamic()
        self.runTest("test_DCB_3D_dyn_coh",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[('\*Variable mass scaling.*', '**'),],   # No mass scaling
                     expectedpyName="test_DCB_3D_dynamic"
                    )

    def test_DCB_3D_dyn_constit_thk(self):
        """ Dynamic DCB, using consitutive thickness != 1 """
        make_dynamic()
        self.runTest("test_DCB_3D_dyn_constit_thk",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[('\*Variable mass scaling.*', '**'),  # No mass scaling
                                    ('constitutive_thk = 1.0', 'constitutive_thk = calc_ct'),
                                    ('density_interface = density', 'density_interface = calc_rho'),
                                    ],
                     expectedpyName="test_DCB_3D_dynamic"
                    )

    def test_DCB_3D_dyn_finite_thk(self):
        """ Dynamic DCB example, using finite thickness cohesive elements """
        make_dynamic()
        interface_el_thk = 0.02
        density_coh = self.param['density']*interface_el_thk
        make_finite_thickness(interface_el_thk, self.param)
        self.runTest("test_DCB_3D_dyn_finite_thk",
                     inpName="test_DCB_3D_compdam",
                     substitutions=[('\*Variable mass scaling.*', '**'),  # No mass scaling
                                    ('interface_el_thk = 0.0', 'interface_el_thk = {}'.format(interface_el_thk)), 
                                    ('density_interface = density', 'density_interface = {}'.format(density_coh)),
                                ],
                     expectedpyName="test_DCB_3D_dynamic",
                     expected_substitutions=[('130000', '306000')]  # Finite thickness requires more increments
                    )

class DCBFatigue(av.TestCase, metaclass=av.ParametricMetaClass):
    """
    Demonstrates the cohesive fatigue model for DCB under displacement control
    """

    # Refers to the template input file name
    baseName = "test_DCB_3D_fatigue"

    # disp is the maximum applied displacement during each fatigue cycle
    parameters = {'disp': [1.48, 1.70, 1.92, 2.25]}

    # staticLoad is the maximum reaction force in the initial static loading step
    expectedpy_parameters = {'staticLoad': [1.813, 2.083, 2.352, 2.754]}

    @classmethod
    def setUpClass(cls):
        av_opts = av.getOptions()
        if av_opts.genPes:
            import unittest
            raise unittest.SkipTest("Skipping pes file generation for DCBFatigue")
        av.copyFiles("_DCB_3D_geometry+sec.inp")

if __name__ == "__main__":
    subroutine_path = os.path.join(os.pardir, os.pardir, 'for', 'CompDam_DGD')
    av.runTests(relPathToUserSub=subroutine_path, double=True)
    
