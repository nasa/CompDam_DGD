parameters = {
	"results": [
		{
            "type": "min",
            "step": "Compression",
            "identifier":
                {
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": -10.0,
            "tolerance": 0.01
        },
		{
            "type": "max",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 78.15,  # ST - 10.0/tan(2.0*alpha0)
            "tolerance": 0.07815
        },
        {
            "type": "finalValue",
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 10.0*0.3,  # Compressive load * coefficient of friction
            "tolerance": 0.03
        },
        {
            "type": "max",
            "step": "Shear",
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
            "step": "Shear",
            "identifier":
                {
                    "symbol": "S23",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.2
        }
	]
}
