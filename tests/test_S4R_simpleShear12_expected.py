parameters = {
	"results": [
		{
            "type": "max",
            "identifier": 
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                },
            "referenceValue": 92.3,
            "tolerance": 0.05
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "Y+",
                    "position": "Node 3"
                },
                { # y
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                }
            ],
            "zeroTol": 0.006,  # Defines how close to zero the y value needs to be
            "referenceValue": 0.0170741,
            "tolerance": 1e-5
        },
        {
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier": 
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        }
	]
}