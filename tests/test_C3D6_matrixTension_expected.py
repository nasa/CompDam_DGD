from utilities import get_enerElas

E2 = 9080.  # Young's modulus, matrix direction
YT = 50.  # Matrix tensile strength
GYT = 0.277  # Mode I fracture toughness
length = 0.1  # Element edge length

enerElas = 0.5*(0.5*YT)**2/E2 * length**3  # computed at stress = 0.5*YT for consistency with where maximum elastic energy in cohesive occurs
enerElas += get_enerElas(YT, GYT, length**2, pen=1e6, dmg_area=0.5)
enerFrac = GYT * length**2

parameters = {
	"results": [
		{
            "type": "max",
            "identifier": 
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": YT,
            "tolerance": YT * 0.001  # 0.1% error
        },
		{
            "type": "max",
            "identifier": 
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                },
            "referenceValue": YT,
            "tolerance": YT * 0.001  # 0.1% error
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U2",
                    "nset": "Y+",
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
            "referenceValue": 2. * GYT / YT,  # Cohesive displacement-jump at complete failure
            "tolerance": 2. * GYT / YT * 0.005  # 0.5% error
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U2",
                    "nset": "Y+",
                    "position": "Node 3"
                },
                { # y
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                }
            ],
            "window": [0.01, 0.015],
            "zeroTol": 0.005,  # Defines how close to zero the y value needs to be
            "referenceValue": 2. * GYT / YT,  # Cohesive displacement-jump at complete failure
            "tolerance": 2. * GYT / YT * 0.005  # 0.5% error
        },
        {
            "type": "max",
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
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                },
            "referenceValue": 1.0,
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
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
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
            "type": "max",
            "identifier": 
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
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
                    "symbol": "SDV_CDM_alpha",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier": 
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "continuous",
            "identifier": 
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 2 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT * area
            "tolerance": enerFrac * 0.002  # 0.02% error
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT * area
            "tolerance": enerFrac * 0.002  # 0.02% error
        },
        {
            "type": "max",
            "identifier": "Strain energy: ALLSE for Whole Model",   # Recoverable strain energy
            "referenceValue": enerElas,  # Elastic strain energy * volume
            "tolerance": enerElas * 0.05
        }
	]
}