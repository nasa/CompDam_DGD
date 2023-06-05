SL = 92.3
GSL = 0.788
length = 0.1

enerFrac = GSL * (length * length)

parameters = {
	"results": [
		{
            "type": "max",
            "identifier": 
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": SL,
            "tolerance": SL * 0.001
        },
		{
            "type": "max",
            "identifier": 
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                },
            "referenceValue": SL,
            "tolerance": SL * 0.001
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
                    "position": "Element 1 Int Point 1"
                }
            ],
            "zeroTol": 0.006,  # Defines how close to zero the y value needs to be
            "referenceValue": 2. * GSL / SL,
            "tolerance": (2. * GSL / SL) * 0.001
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
                    "position": "Element 2 Int Point 1"
                }
            ],
            "zeroTol": 0.006,  # Defines how close to zero the y value needs to be
            "referenceValue": 2. * GSL / SL,
            "tolerance": (2. * GSL / SL) * 0.001
        },
        {
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
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
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
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
                    "position": "Element 1 Int Point 1"
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
                    "position": "Element 2 Int Point 1"
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
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "continuous",
            "identifier": 
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GSL*LC1*LC3
            "tolerance": enerFrac * 0.001
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GSL*LC1*LC3
            "tolerance": enerFrac * 0.001
        }
	]
}