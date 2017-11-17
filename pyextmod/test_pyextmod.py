import numpy as np
import CompDam_DGD
import viz

# logging
CompDam_DGD.dgd_mod.log_init(level=4, filename='pyextmod_run_output.txt')

# Load material properties from props file
m = CompDam_DGD.matprop_mod.loadmatprops('IM7-8552', False, 4, [10300, 0, 0.2, 0])

# Load parameters from the 'CompDam.parameters' file
p = CompDam_DGD.parameters_mod.loadparameters()

# State variables
svarray = [
	0,
	0.243194794854592E+00,
	-.935227656912301E+00,
	0.753295466440013E-03,
	0,
	0.199999800000001E+00,
	0.200000050000001E+00,
	0.182999968000000E+00,
	0,
	0,
	0.140129846432482E-43,
	-.910526921593782E+00,
	0.910527510918844E+00,
	0.100000000000000E+01,
	0.000000000000000E+00,
	0.000000000000000E+00,
	0.100000000000000E+01,
	0,
	0.100000000000000E-05,
	-.362036217196037E-01,
	0.583479715694666E+00
]
sv = CompDam_DGD.statevar_mod.loadstatevars(len(svarray), svarray, m)
# print sv


# Other inputs
F = np.array([
	[0.825492123341086E+00, 0.757572103784855E-03, 0.301104249035103E-03],
	[-.146788494201401E+00, 0.100523463166065E+01, 0.641873185909126E-03],
	[-.662872667983369E-03, 0.871848067078203E-04, 0.100738051937547E+01]])

F_old = np.array([
	[0.828052171305039E+00, 0.741124652528062E-03, 0.313281227086903E-03],
	[-.148260374255618E+00, 0.100527363842522E+01, 0.620850904202262E-03],
	[-.624763825230815E-03, 0.790313763395604E-04, 0.100728719006083E+01]])

U = np.array([
	[0.834616447357706E+00, -.799991473742500E-01, -.263339433701766E-03],
	[-.799991473742500E-01, 0.100204653434751E+01, 0.354440612331553E-03],
	[-.263339433701766E-03, 0.354440612331553E-03, 0.100738067209290E+01]])

Cauchy = np.zeros((3,3), order='F')

CompDam_DGD.dgd_mod.dgdevolve(u=U, f=F, f_old=F_old, m=m, p=p, sv=sv, ndir=3, nshr=3, dt=0, cauchy=Cauchy)
# print Cauchy

CompDam_DGD.dgd_mod.log_close()


# Visualize
viz.visualize(lc=svarray[5:8], alpha=svarray[9], F=F, pathToLogFile='pyextmod_run_output.txt', initMD=2, initEQk=1)
