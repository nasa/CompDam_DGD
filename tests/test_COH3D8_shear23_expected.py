import math

SL = 92.3  # longitudinal shear strength
GSL = 0.788  # Mode II matrix fracture toughness
YC = 199.8  # matrix compressive strength
alpha0 = 0.925  # matrix crack orientation due to pure matrix compression failure
length = 0.2  # element edge length

ST = YC * math.cos(alpha0) * (math.sin(alpha0) + math.cos(alpha0) / math.tan(2*alpha0))
enerFrac = GSL * length * length

parameters = {
	"results": [
		{
            "type": "max",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": ST,
            "tolerance": ST * 0.001  # 0.1% error
        },
        {
            "type": "disp_at_zero_y",
            "step": "Shear",
            "identifier": [
                { # x
                    "symbol": "U2",
                    "nset": "Z+",
                    "position": "Node 5"
                },
                { # y
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.015, 0.025],
            "zeroTol": 0.005,  # Defines how close to zero the y value needs to be
            "referenceValue": 2.0 * GSL / ST,  # Cohesive displacement-jump at complete failure
            "tolerance": 2.0 * GSL / ST * 0.001  # 0.1% error
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
            "tolerance": 0.1
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
