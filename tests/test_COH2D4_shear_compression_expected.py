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
            "referenceValue": -10.0,
            "tolerance": 0.01
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
            "referenceValue": 95.8,  # SL - 10.0 * SL/YC*cos(2.0*alpha0)/cos(alpha0)**2
            "tolerance": 0.0958
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
            "referenceValue": 2.0*0.788/95.8,  # u_f = 2*GSL/(SL - 10.0 * SL/YC*cos(2.0*alpha0)/cos(alpha0)**2)
            "tolerance": 1e-5
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
        }
	]
}
