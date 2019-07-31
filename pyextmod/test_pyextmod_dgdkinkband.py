import numpy as np
import CompDam_DGD

# logging
CompDam_DGD.dgd_mod.log_init(level=4, filename='pyextmod_run_output.txt')

# Load material properties from props file
# m = CompDam_DGD.matprop_mod.loadmatprops('IM7-8552', False, 4, [10300, 0, 0.2, 0])
m = CompDam_DGD.matprop_mod.loadmatprops('IM7-8552', 40, 
	[10300,        0,       0.2,       0,        0,       0,          0,      0, 
	140653.,    8703.,    5164.,     0.32,     0.45,     80.1,     97.6,     0.24,
	0.739,     2.07,    288.2,    0.925,    8703.,    5164.,    3001.,     0.32,
	-5.5e-06, 2.58e-05, 4.06e-09,      5.4,   2326.2,      0.2,     205.,      0.5,
	1730.6,      0.2,      61.,      0.5,      10.,     0.05,       0.,      0.3])
CompDam_DGD.matprop_mod.consistencychecks(m, True)


# Load parameters from the 'CompDam.parameters' file
p = CompDam_DGD.parameters_mod.loadparameters()

# State variables
svarray = [
	0,                        # d2
	0.100010907174600E+01,    # Fb1
	0.285074113899452E-02,    # Fb2
	-.564065053387765E-03,    # Fb3
	0,                        # B
	0.251999888476518E+00,    # LC1
	0.252000390522334E+00,    # LC2
	0.182999998000000E+00,    # LC3
	0,                        # FIm
	0,                        # alpha
	1,                        # STATUS
	0.834099516228293E-03,    # Plas12
	0.834106365248572E-03,    # Inel12
	0.100000000000000E+01,    # FIft
	0.000000000000000E+00,    # slide 1
	0.000000000000000E+00,    # slide 2
	0.393443624534415E-02,    # FIfc
	0,                        # d1T
	0.100000000000000E-05,    # d1C
	0.000000000000000E+00,    # Plas13
	0.000000000000000E+00,    # Inel13
	-.121216129159982E-01,    # phi0
	0.288105568001332E-02,    # gamma
	0.000000000000000E+00,    # Fm1
	0.000000000000000E+00,    # Fm2
	0.000000000000000E+00     # Fm3
]
sv = CompDam_DGD.statevar_mod.loadstatevars(len(svarray), svarray, m)
print(sv)


# Other inputs
F = np.array([
	[0.100010537052244E+01, 0.671448581859227E-02, -.394071455445696E-04],
	[0.215218597548654E-02, 0.995441501182020E+00, -.100155558582461E-02],
	[-.559184398241147E-03, -.122772703106108E-02, 0.100136989172593E+01]])

F_old = np.array([
	[0.100010292495566E+01, 0.669475133238129E-02, -.355545830703892E-04],
	[0.217728430299811E-02, 0.995457191114569E+00, -.995635844788960E-03],
	[-.564742609464082E-03, -.122174252924641E-02, 0.100137225125743E+01]])

U = np.array([
	[0.100009794753325E+01, 0.443883268741269E-02, -.298066487405907E-03],
	[0.443883268741269E-02, 0.995454383031751E+00, -.111443987390743E-02],
	[-.298066487405907E-03, -.111443987390743E-02, 0.100136972887234E+01]])

Cauchy = np.zeros((3,3), order='F')

enerintern = 0
enerinelas = 0

CompDam_DGD.dgd_mod.dgdkinkband(u=U, f=F, f_old=F_old, m=m, p=p, sv=sv, ndir=3, nshr=3, dt=0, density_abq=1, cauchy=Cauchy, enerintern=enerintern, enerinelas=enerinelas)
# print Cauchy

CompDam_DGD.dgd_mod.log_close()
