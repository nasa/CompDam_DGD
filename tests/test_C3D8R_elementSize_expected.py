Lc1 = 0
Lc2 = 0

parameters = {
	"results": [
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_Lc2",
                    "elset": "SINGLEELEMENT",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": Lc2,
            "tolerance": Lc2*0.001
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_Lc3",
                    "elset": "SINGLEELEMENT",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.1,
            "tolerance": 0.0001
        }
	]
}
