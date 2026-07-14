from utilities import etaL, get_enerElas

applied_compression = -10.0

SL = 92.3  # longitudinal shear strength
YC = 199.8  # matrix compressive strength
GSL = 0.788  # Mode II matrix fracture toughness
alpha0 = 0.925  # matrix crack orientation due to pure matrix compression failure
length = 0.2  # element edge length
out_of_plane_thk = 0.5

eta_L = etaL(SL, YC, alpha0)
enerElas_dmg_0p5 = get_enerElas(SL, GSL, length*out_of_plane_thk)

enerFrac = GSL * length * out_of_plane_thk

parameters = {
	"results": [
		{
            "type": "min",
            "step": "Compression",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": applied_compression,
            "tolerance": applied_compression * 0.001  # 0.1% error
        },
		{
            "type": "max",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": SL - eta_L * applied_compression,
            "tolerance": (SL - eta_L * applied_compression) * 0.001  # 0.1% error
        },
        {
            "type": "disp_at_zero_y",
            "step": "Shear",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "Y+",
                    "position": "Node 4"
                },
                { # y
                    "symbol": "S12",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.015, 0.020],
            "zeroTol": 0.005,  # Defines how close to zero the y value needs to be
            "referenceValue": 2. * GSL / (SL - eta_L * applied_compression),
            "tolerance": 2. * GSL / (SL - eta_L * applied_compression) * 0.001  # 0.1% error
        },
        {
            "type": "max",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "SDV_COH_dmg",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.2
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GSL*area
            "tolerance": enerFrac * 0.001  # 0.1% error
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GSL*area
            "tolerance": enerFrac * 0.001  # 0.1% error
        },
        {
            "type": "max",
            "identifier": "Strain energy: ALLSE for Whole Model",  # Recoverable strain energy
            "referenceValue": enerElas_dmg_0p5,  # Elastic strain energy * volume
            "tolerance": enerElas_dmg_0p5 * 0.002
        }
	]
}
