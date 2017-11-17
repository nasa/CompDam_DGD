parameters = {
	"results": [
		{
			"type": "log_stress_at_failure_init",
			"step": "Step-1",
			"failureIndices": [
				{
					"symbol": "SDV_CDM_FIM",
					"elset": "ALL",
					"position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
				}
			],
			"stressComponents": [
				{
					"symbol": "S12",
					"elset": "ALL",
					"position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
				},
				{
					"symbol": "S22",
					"elset": "ALL",
					"position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
				}
			],
			"additionalIdentifiersToStore": [
				{
					"symbol": "SDV_CDM_ALPHA",
					"elset": "ALL",
					"position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
				}
			]
		},
		{
			"type": "max",
			"step": "Step-1",
			"identifier": {
                "symbol": "SDV_CDM_FIM",
                "elset": "ALL",
                "position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
            },
			"referenceValue": 1.0,
			"tolerance": 0.0
		}
	]
}
