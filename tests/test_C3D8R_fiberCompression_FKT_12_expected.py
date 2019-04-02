# Variable results
crushStress=-9.55

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
            "referenceValue": -36.003,
            "tolerance": 0.72  # 2%
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
            "referenceValue": -0.00157,
            "tolerance": 0.0001
        }
	]
}
