stress_ratio = 0.5
strength = 50.
toughness = 0.277
modulus = 9080.0
length = 0.1

delta_i = strength / (modulus / length)
delta_f = 2.0 * toughness / strength
delta_inflection_fatigue = delta_f - stress_ratio * (delta_f - delta_i)

parameters = {
	"results": [
		{
            "type": "max",
			"step": "Load",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": strength * stress_ratio,
            "tolerance": strength * stress_ratio * 0.002
        },
		{
            "type": "max",
			"step": "Fatigue",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": strength * stress_ratio,
            "tolerance": strength * stress_ratio * 0.003
        },
        {
            "type": "finalValue",
            "identifier": 
                {
                    "symbol": "SDV_CDM_alpha",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "xy_infl_pt_bilinear",
            "step": "Fatigue",
            "identifier": [
                { # x
                    "symbol": "U2",
                    "nset": "LOADNODE",
                    "position": "Node 3"
                },
                { # y
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
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
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "disp_at_zero_y",
            "step": "Fatigue",
            "identifier": [
                { # x
                    "symbol": "U2",
                    "nset": "LOADNODE",
                    "position": "Node 3"
                },
                { # y
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.01, 0.015],
            "zeroTol": 0.005,  # Defines how close to zero the y value needs to be
            "referenceValue": 2. * toughness / strength,  # Cohesive displacement-jump at complete failure
            "tolerance": 2. * toughness / strength * 0.05  # 5% error
        },
	]
}
