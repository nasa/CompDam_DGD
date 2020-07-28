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
            "referenceValue": 62.3, # YT
            "tolerance": 0.0623
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
            "referenceValue": 0.00889, # u_f = 2*GYT/YT
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
                    "symbol": "S22",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        }
	]
}
