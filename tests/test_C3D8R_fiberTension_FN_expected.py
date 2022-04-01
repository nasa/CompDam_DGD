parameters = {
    "results": [
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 2326.2,
            "tolerance": 60.0
        },
        {
            "type": "xy_infl_pt",
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
            "window": [0.05, 0.06],
            "referenceValue": [0.05730, 465.24],
            "tolerance": [0.003, 65]
        },
        {
            "type": "disp_at_zero_y",
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
            "window": [0.15, 0.4],  # [min, max]
            "zeroTol": 0.00623,  # Defines how close to zero the y value needs to be
            "referenceValue": 0.2865,
            "tolerance": 0.055
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.2
        }
    ]
}
