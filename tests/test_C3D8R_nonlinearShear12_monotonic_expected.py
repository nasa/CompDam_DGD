tau = 90  # Stress at 5% strain from RO law
aPL = 4.412E-10
nPL = 5.934
G12 = 5290.
enerElas = 0.5*tau**2/G12*0.1*0.1*0.1
enerPlas = aPL/G12*tau**(nPL + 1.)*nPL/(nPL + 1.)*0.1**3

parameters = {
	"results": [
		{
            "type": "max",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": tau,
            "tolerance": 0.01*tau
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "max",
            "identifier": "Strain energy: ALLSE for Whole Model",   # Recoverable strain energy
            "referenceValue": enerElas,   # Elastic strain energy * volume: 0.5*stress^2/G12*LC1*LC2*LC3
            "tolerance": 0.02*enerElas
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerPlas,  # Plastic energy dissipation from RO law * volume: aPL/Modulus*tau**(nPL + one)*nPL/(nPL + one)*LC1*LC2*LC3
            "tolerance": 0.02*enerPlas  # Enfore 2% error maximum. Note most of the error is related to large strain.
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerPlas,  # Plastic energy dissipation from RO law * volume: aPL/Modulus*tau**(nPL + one)*nPL/(nPL + one)*LC1*LC2*LC3
            "tolerance": 0.02*enerPlas
        }
	]
}