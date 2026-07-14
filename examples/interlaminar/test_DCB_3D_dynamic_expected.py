peakLoad = 8.0
finalLoad = 4.38

parameters = {
	"results": [
        {
            "type": "tabular",
            "identifier": {
                "symbol": "RF3",
                "nset": "T_NODE"
            },
            "referenceValue": [
                            (2.0e-4, 0.031),
                            (4.2e-4, 0.1),
                            (5.2e-4, 0.17),
                            ],
            "tolerance_percentage": 0.5
        },
        {
            "type": "length",
            "identifier": {
                "symbol": "RF3",
                "nset": "T_NODE"
            },
            "referenceValue": 130000,
            "tolerance": 10000
        },
	]
}
