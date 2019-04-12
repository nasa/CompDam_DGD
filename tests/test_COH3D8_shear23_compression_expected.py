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
            "referenceValue": -10.0,
            "tolerance": 0.01
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
            "referenceValue": 78.15,  # ST - 10.0/tan(2.0*alpha0)
            "tolerance": 0.07815
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
            "referenceValue": 2.0*0.788/78.15,  # u_f = 2*GSL/(ST - 10.0/tan(2.0*alpha0))
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
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        }
	]
}
