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
            "referenceValue": -71,
            "tolerance": 1.0
        },
        {
            "type": "max",
            "step": "Step-1",
            "identifier": {
                "symbol": "SDV_CDM_PLAS12",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.137,
            "tolerance": 0.01
        }
    ]
}
