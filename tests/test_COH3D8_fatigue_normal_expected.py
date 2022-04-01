stress_ratio = 0.5
strength = 80.1
toughness = 0.240
penalty = 2.e5

delta_i = strength / penalty
delta_f = 2.0 * toughness / strength
delta_inflection_fatigue = delta_f - stress_ratio * (delta_f - delta_i)

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
            "type": "xy_infl_pt_bilinear",
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
            "window": [
				delta_inflection_fatigue - (delta_inflection_fatigue - delta_i) * 0.80,
				delta_inflection_fatigue + (delta_f - delta_inflection_fatigue) * 0.80
				],
            "referenceValue": [delta_inflection_fatigue, strength * stress_ratio],
            "tolerance": [delta_inflection_fatigue * 0.02, strength * stress_ratio * 0.02],
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
        }
	]
}
