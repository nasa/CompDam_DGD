import numpy as np
import CompDam_DGD
import helpers as h

CompDam_DGD.dgd_mod.log_init(level=4, filename='pyextmod_run_output.txt')
(m, p, sv, sv_old, F, F_old, U, debugpy) = h.loaddebugpy(filename='<debug-file-name-nodotpy>')

# Run CompDam
sv_calculated  = sv_old
Cauchy = np.zeros((3,3), order='F')
func = getattr(CompDam_DGD.dgd_mod, debugpy.called_from.lower())
func(u=U, f=F, f_old=F_old, m=m, p=p, sv=sv_calculated, ndir=3, nshr=3, dt=0, cauchy=Cauchy)

CompDam_DGD.dgd_mod.log_close()
