parameters = {
	"results": [
		{
            "type": "max",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 92.3,  # SL
            "tolerance": 0.0923
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
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
            "referenceValue": 2.0*0.788/92.3,  # u_f = 2*GSL/SL
            "tolerance": 1e-5
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
                    "symbol": "S12",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.2
        }
	]
}
