parameters = {
	"results": [
        {
            "type": "slope",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "LE22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
                { # y
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.0001, 0.005],  # [min, max] in x
            "referenceValue": 9080, # E2
            "tolerance": 90.8 # 1%
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_FIm",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        }
	]
}
