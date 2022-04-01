class Materials():
	'''
	Published material property sets for verification and demonstration
	'''

	materials = {}

	# Original from Soden et al 1998 (WWFE-I). Not sure why some values in LaRC03 and LaRC04 differ from Soden.
	EGlassLY556_larc03 = {
		'E1': 53480.,
		'E2': 17700.,
		'nu12': 0.3,
		'G12': 5830.,
		'alpha_PL': 8.22e-13,
		'eta_PL': 7.66,
		'XT': 1140.,
		'XC': 570.,
		'SL': 61.,
		'YC': 138.,
		'YT': 39,   # This is the number I get from the graph; not sure why the table has 36
		'alpha0': 0.925
	}
	materials['EGlassLY556_larc03'] = EGlassLY556_larc03
	EGlassLY556_larc04 = {
		'E1': 53480.,
		'E2': 17700.,
		'nu12': 0.3,
		'G12': 5830.,
		'alpha_PL': 8.22e-13,
		'eta_PL': 7.66,
		'XT': 1140.,
		'XC': 570.,
		'SL': 66.5,
		'YC': 130.3,
		'YT': 37.5,
		'alpha0': 0.925
	}
	materials['EGlassLY556_larc04'] = EGlassLY556_larc04

	EGlassMY750 = {
		'E1': 45600.,
		'E2': 16200.,
		'nu12': 0.278,
		'G12': 5830.,
		'alpha_PL': 3.253e-13,
		'eta_PL': 7.88,
		'XT': 1280.,
		'XC': 800.,
		'SL': 90.,
		# 'SL': 73.,
		'YC': 145.,
		'YT': 40.,
		'alpha0': 0.925
	}
	materials['EGlassMY750'] = EGlassMY750

	IM78552 = {
		'E1': 171420.,
		'E2': 9080.,
		'nu12': 0.32,
		'G12': 5290.,
		'alpha_PL': 4.412e-10,
		'eta_PL': 5.934,
		'XT': 2326.2,
		'XC': 1200.1,
		'SL': 92.3,
		'YC': 199.8,
		'YT': 62.3,
		'GYT': 0.2770,
		'GSL': 0.7880,
		'eta_BK': 1.634,
		'alpha0': 0.925
	}
	materials['IM78552'] = IM78552

	IM78552ACC = {
		'E1': 140653.,
		'E2': 8703.,
		'nu12': 0.32,
		'G12': 5164.,
		'alpha_PL': 4.06e-09,
		'eta_PL': 5.4,
		'XT': 2326.2,
		'XC': 1730.6,
		'SL': 97.6,
		'YC': 288.2,
		'YT': 80.1,
		'GYT': 0.240,
		'GSL': 0.7390,
		'eta_BK': 2.07,
		'alpha0': 0.925
	}
	materials['IM78552ACC'] = IM78552ACC
