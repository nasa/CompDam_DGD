loadRatio = 0.0
matrixStrength = 62.3
fiberStrength = 2326.2
s11_expected = fiberStrength if loadRatio > 0.5 else 2*loadRatio*fiberStrength
s22_expected = matrixStrength if loadRatio <= 0.5 else 2*(1-loadRatio)*matrixStrength

parameters = {
	"results": [
        {
            "type": "finalValue",
			"step": "Step-1",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": s11_expected,
            "tolerance": abs(fiberStrength * 0.05)
        },
        {
            "type": "finalValue",
			"step": "Step-1",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": s22_expected,
            "tolerance": abs(matrixStrength * 0.05)
        },
		{
			"type": "log_stress_at_failure_init",
			"step": "Step-1",
			"failureIndices": [
				{
					"symbol": "SDV_CDM_FIM",
					"elset": "ALL",
					"position": "Element 1 Int Point 1"
				},
				{
					"symbol": "SDV_CDM_FIFT",
					"elset": "ALL",
					"position": "Element 1 Int Point 1"
				},
				{
					"symbol": "SDV_CDM_FIFC",
					"elset": "ALL",
					"position": "Element 1 Int Point 1"
				}
			],
			"stressComponents": [
				{
					"symbol": "S11",
					"elset": "ALL",
					"position": "Element 1 Int Point 1"
				},
				{
					"symbol": "S22",
					"elset": "ALL",
					"position": "Element 1 Int Point 1"
				}
			],
			"additionalIdentifiersToStore": [
				{
					"symbol": "SDV_CDM_ALPHA",
					"elset": "ALL",
					"position": "Element 1 Int Point 1"
				}
			]
		}
	]
}
