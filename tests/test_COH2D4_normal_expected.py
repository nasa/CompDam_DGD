YT = 62.3
GYT = 0.277
length = 0.2

enerFrac = GYT * length * 1.0

parameters = {
	"results": [
		{
            "type": "max",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": YT,
            "tolerance": YT * 0.001  # 0.1% error
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U2",
                    "nset": "Y+",
                    "position": "Node 4"
                },
                { # y
                    "symbol": "S22",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.008, 0.010],
            "zeroTol": 0.005,  # Defines how close to zero the y value needs to be
            "referenceValue": 2. * GYT / YT,  # Cohesive displacement-jump at complete failure
            "tolerance": 2. * GYT / YT * 0.001  # 0.1% error
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
                    "symbol": "S22",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT*area
            "tolerance": enerFrac * 0.001  # 0.1% error
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT*area
            "tolerance": enerFrac * 0.001  # 0.1% error
        }
	]
}
