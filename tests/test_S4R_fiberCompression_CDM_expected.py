parameters = {
	"results": [
        {
            "type": "min",
            "step": "Step-1",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
            },
            "referenceValue": -1200.1,
            "tolerance": 20.0
        },
        {
            "type": "xy_infl_pt",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                }
            ],
            "window": [-0.06, -0.03],
            "referenceValue": [-0.03958, -240.02],
            "tolerance": [0.0025, 40]
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                }
            ],
            "window": [-0.35, -0.1],  # [min, max]
            "zeroTol": 0.00623,  # Defines how close to zero the y value needs to be
            "referenceValue": -0.1979,
            "tolerance": 0.075
        },
        {
            "type": "max",
            "step": "Step-1",
            "identifier": {
                "symbol": "SDV_CDM_D2",
                "elset": "ALL",
                "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
            },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
            },
            "referenceValue": 0.0,
            "tolerance": 0.9
        }
	]
}
