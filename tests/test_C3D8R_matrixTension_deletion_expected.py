from utilities import get_enerElas

E2 = 9080.  # Young's modulus, matrix direction
YT = 62.3   # Matrix tensile strength
GYT = 0.277  # Mode I fracture toughness
lengthx = 0.1
lengthy = 0.15
lengthz = 0.2

disp_at_zt = 2. * GYT / YT  # Cohesive displacement-jump at complete failure

enerElas = 0.5*(0.5*YT)**2/E2 * lengthx*lengthy*lengthz  # computed at stress = 0.5*YT for consistency with where maximum elastic energy in cohesive occurs
enerElas += get_enerElas(YT, GYT, lengthx*lengthz, pen=1e6, dmg_area=0.5)
enerFrac = GYT * lengthx*lengthz

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
            "window": [0.8*disp_at_zt, 1.2*disp_at_zt],
            "zeroTol": 0.005,  # Defines how close to zero the y value needs to be
            "referenceValue": disp_at_zt,
            "tolerance": 2. * GYT / YT * 0.001  # 0.1% error
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
            "type": "finalValue",
            "identifier":
                {
                    "symbol": "SDV_CDM_STATUS",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {  # This is a proxy for element deletion; after the element is deleted density does not change
            "type": "finalValue",
            "identifier": "Density: DENSITY at Element 1 Int Point 1 in ELSET ALL_ELEMS",
            "referenceValue": 1.48e-5,
            "tolerance": 1.48e-5 * 0.01
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
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT*LC1*LC3
            "tolerance": enerFrac * 0.002  # 0.2% error
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT*LC1*LC3
            "tolerance": enerFrac * 0.002  # 0.2% error
        },
        {
            "type": "max",
            "identifier": "Strain energy: ALLSE for Whole Model",  # Recoverable strain energy
            "referenceValue": enerElas,  # Elastic strain energy * volume
            "tolerance": enerElas * 0.04
        }
	]
}