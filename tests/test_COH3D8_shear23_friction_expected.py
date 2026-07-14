from utilities import etaT, get_ST, get_enerElas

applied_compression = -10.0
coefficient_of_friction = 0.3

SL = 92.3  # longitudinal shear strength
YC = 199.8  # matrix compressive strength
GSL = 0.788  # Mode II matrix fracture toughness
alpha0 = 0.925  # matrix crack orientation due to pure matrix compression failure
length = 0.2  # element edge length

ST = get_ST(YC, alpha0)
eta_T = etaT(alpha0)
enerElas_dmg_0p5 = get_enerElas(SL, GSL, length**2)

enerFrac = GSL * length * length

parameters = {
	"results": [
		{
            "type": "min",
            "step": "Compression",
            "identifier":
                {
                    "symbol": "S33",
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
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": ST - eta_T * applied_compression,
            "tolerance": (ST - eta_T * applied_compression) * 0.001  # 0.1% error
        },
        {
            "type": "finalValue",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": -applied_compression * coefficient_of_friction,
            "tolerance": -applied_compression * coefficient_of_friction * 0.001  # 0.1% error
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
                    "symbol": "S23",
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
            "tolerance": enerElas_dmg_0p5 * 0.05  # 5% error here allowed because the reference value does not account for friction
        }
	]
}
