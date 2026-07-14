import math

def get_enerElas(strength, toughness, area, pen=1e6, dmg_area=0.5):
    disp_initiation = strength/pen
    disp_final = 2*toughness/strength
    dmg_pen = dmg_area/(dmg_area + (disp_initiation/disp_final)*(1 - dmg_area)) # Penalty damage variable when area damage variable = 0.5
    disp = disp_initiation+(disp_final-disp_initiation)*dmg_area
    stress = pen*(1-dmg_pen)*disp
    return 0.5*disp*stress * area

def etaL(SL, YC, alpha0):
    return -SL * math.cos(2*alpha0) / (YC * math.cos(alpha0) * math.cos(alpha0))

def etaT(alpha0):
    return -1./math.tan(2*alpha0)

def get_ST(YC, alpha0):
    return YC * math.cos(alpha0) * (math.sin(alpha0) + math.cos(alpha0) / math.tan(2*alpha0))
