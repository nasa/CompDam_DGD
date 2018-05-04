# Variable results
crushStress=-13.1

parameters = {
	"results": [
        {
            "type": "min",
            "step": "Step-1",
            "identifier": {
                "symbol": "RF1",
                "nset": "LOADAPP",
                "position": "Node 4"
            },
            "referenceValue": -29.3,
            "tolerance": 0.586  # 2%
        },
        {
            "type": "finalValue",
            "step": "Step-1",
            "identifier": {
                "symbol": "RF1",
                "nset": "LOADAPP",
                "position": "Node 4"
            },
            "referenceValue": crushStress,
            "tolerance": 0.02*abs(crushStress)  # 2%
        },
        {
            "type": "x_at_peak_in_xy",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "RF1",
                    "nset": "LOADAPP"
                }
            ],
            "referenceValue": -0.00141,
            "tolerance": 0.00001
        }
	]
}
