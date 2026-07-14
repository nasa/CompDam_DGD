peakLoad = 8.0
peakLoadTol = 0.03
finalLoad = 4.38
finalLoadTol = 0.05
step_time = 1.0
freq = 200/step_time

parameters = {
	"results": [
        {
            "type": "max",
			"step": "Load",
            "identifier": {
                "symbol": "RF3",
                "nset": "T_NODE",
                "position": "Node 9999998",
                "filterCutOffFreq": freq
            },
            "referenceValue": peakLoad,
            "tolerance": peakLoad * peakLoadTol,
        },
        {
            "type": "finalValue",
            "identifier": {
                "symbol": "RF3",
                "nset": "T_NODE",
                "position": "Node 9999998",
                "filterCutOffFreq": freq
            },
            "referenceValue": finalLoad,
            "tolerance": finalLoad * finalLoadTol,
        },
	]
}
