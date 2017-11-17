parameters = {
	"results": [
        {
            "type": "continuous",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.0,
            "tolerance": 0.5
        },
        {
            "type": "max",
            "step": "Step-1",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 2326.2,
            "tolerance": 60.0
        },
        {
            "type": "min",
            "step": "Step-3",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": -1200.1,
            "tolerance": 20.0
        },
        {
            "type": "max",
            "step": "Step-5",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 1273,
            "tolerance": 5.0
        },
        {
            "type": "min",
            "step": "Step-7",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": -745,
            "tolerance": 5.0
        },
        {
            "type": "max",
            "step": "Step-9",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 410,
            "tolerance": 5.0
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-9",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.15, 0.45],  # [min, max]
            "zeroTol": 0.00623,  # Defines how close to zero the y value needs to be
            "referenceValue": 0.2865,
            "tolerance": 0.055
        },
        {
            "type": "max",
            "step": "Step-9",
            "identifier": {
                "symbol": "SDV_CDM_D2",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.0,
            "tolerance": 0.0
        }
	]
}
