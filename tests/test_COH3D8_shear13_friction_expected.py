import math

applied_compression = -10.0
coefficient_of_friction = 0.3

SL = 92.3  # longitudinal shear strength
YC = 199.8  # matrix compressive strength
GSL = 0.788  # Mode II matrix fracture toughness
alpha0 = 0.925  # matrix crack orientation due to pure matrix compression failure
length = 0.2  # element edge length

eta_L = -SL * math.cos(2*alpha0) / (YC * math.cos(alpha0) * math.cos(alpha0))

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
                    "symbol": "S13",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": SL - eta_L * applied_compression,
            "tolerance": (SL - eta_L * applied_compression) * 0.001  # 0.1% error
        },
        {
            "type": "finalValue",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S13",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": -applied_compression * coefficient_of_friction,
            "tolerance": -applied_compression * coefficient_of_friction * 0.001  # 0.1%
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
                    "symbol": "S13",
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
        }
	]
}
