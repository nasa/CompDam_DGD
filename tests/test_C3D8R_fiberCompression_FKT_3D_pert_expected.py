rf1 = -1*1200.1*(0.15*0.15)
crushStress = -7.9
dispatpeakload = -0.001326

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
            "referenceValue": rf1,
            "tolerance": 0.02*abs(rf1)  # 2%
        },
        {
            "type": "finalValue",
            "step": "Step-1",
            "identifier": {
                "symbol": "SDV_CDM_GAMMA_13",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.325,
            "tolerance": 0.003
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
            "referenceValue": dispatpeakload,
            "tolerance": 0.0001
        }
	]
}
