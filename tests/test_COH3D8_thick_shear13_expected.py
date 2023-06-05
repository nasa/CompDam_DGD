SL = 92.3
GSL = 0.788
length = 0.2

enerFrac = GSL * length * length

parameters = {
	"results": [
		{
            "type": "max",
            "identifier":
                {
                    "symbol": "S13",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": SL,
            "tolerance": SL * 0.001  # 0.1% error
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "Z+",
                    "position": "Node 5"
                },
                { # y
                    "symbol": "S13",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.015, 0.020],
            "zeroTol": 0.005,  # Defines how close to zero the y value needs to be
            "referenceValue": 2. * GSL / SL,  # Cohesive displacement-jump at complete failure
            "tolerance": 2. * GSL / SL * 0.001  # 0.1% error
        },
        {
            "type": "max",
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
