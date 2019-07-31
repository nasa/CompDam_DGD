stress_ratio = 0.5
strength = 62.3
toughness = 0.277

delta_f = 2.0 * toughness / strength

parameters = {
	"results": [
		{
            "type": "max",
			"step": "Load",
            "identifier":
                {
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": strength * stress_ratio,
            "tolerance": strength * stress_ratio * 0.001
        },
		{
            "type": "max",
			"step": "Fatigue",
            "identifier":
                {
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": strength * stress_ratio,
            "tolerance": strength * stress_ratio * 0.003
        },
        {
            "type": "disp_at_zero_y",
            "step": "Fatigue",
            "identifier": [
                { # x
                    "symbol": "U3",
                    "nset": "Z+",
                    "position": "Node 5"
                },
                { # y
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [delta_f * 0.80, delta_f * 1.20],
            "zeroTol": strength * 0.001,  # Defines how close to zero the y value needs to be
            "referenceValue": delta_f,
            "tolerance": delta_f * 0.1
        },
        {
            "type": "max",
			"step": "Fatigue",
            "identifier":
                {
                    "symbol": "SDV_COH_dmg",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
			"step": "Fatigue",
            "identifier":
                {
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": strength * 0.1
        }
	]
}
