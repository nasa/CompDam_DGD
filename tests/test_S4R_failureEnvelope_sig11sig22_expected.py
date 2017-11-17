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
				},
				{
					"symbol": "SDV_CDM_FIFT",
					"elset": "ALL",
					"position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
				},
				{
					"symbol": "SDV_CDM_FIFC",
					"elset": "ALL",
					"position": "Element 1 Int Point 1 Sec Pt SPOS, (fraction = 1:0)"
				}
			],
			"stressComponents": [
				{
					"symbol": "S11",
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
		}
	]
}
