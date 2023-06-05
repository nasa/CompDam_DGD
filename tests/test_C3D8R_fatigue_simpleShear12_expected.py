stress_ratio = 0.5
strength = 92.3
toughness = 0.788
modulus = 5290.0
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
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": strength * stress_ratio,
            "tolerance": strength * stress_ratio * 0.005
        },
		{
            "type": "max",
			"step": "Fatigue",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": strength * stress_ratio,
            "tolerance": strength * stress_ratio * 0.005
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
            "type": "finalValue",
            "identifier": 
                {
                    "symbol": "U1",
                    "nset": "LOADNODE",
                    "position": "Node 3"
                },
            "referenceValue": delta_f,
            "tolerance": delta_f * 0.1
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
	]
}