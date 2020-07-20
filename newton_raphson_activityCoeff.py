import gSAFTmm
from scipy import optimize
import math

compounds = [
    "water",
    "Na+",
    "Cl-"
]

T = 298.15

P = 101325

z_ref = [1, 0, 0]

gfNa = -261.905
gfCl = -131.228
gfNaCl = -384.138

################# Define System #########################

compounds_toString = ','.join(compounds)
system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)

################### Ksp getter ##########################

delta_g = (gfNa + gfCl - gfNaCl) * 1000
R = 8.3144626
exponent = - delta_g / (R * T)
Ksp = math.exp(exponent)

print('Ksp: ' + str(Ksp))

################ activity obj func #######################

def activity_product(m):
    prediction = S.LiquidActivityCoefficientAsymMolal(T,P,[55.5093,m,m],z_ref)
    coeff_arr = prediction.values()
    product_of_activity = 1
    for i in range(len(coeff_arr)):
        if i==0: continue
        product_of_activity *= (coeff_arr[i] * m)
    obj_func = Ksp - product_of_activity
    return obj_func

################ newton raphson ###########################

result = optimize.newton(activity_product, 6.88)

print('result: ' + str(result))

