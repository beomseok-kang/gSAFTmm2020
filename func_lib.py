import gSAFTmm
from scipy import optimize
import math

################# Define System #########################

def saft_system(compounds_arr):
    compounds_toString = ','.join(compounds_arr)
    system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"
    return gSAFTmm.System(system_string)

############## H and S getter ###########################

def S_ref_getter(delta_H_ref_val, delta_G_ref_val, temp_ref):
    return (delta_H_ref_val - delta_G_ref_val) / temp_ref

def delta_H(temperature, temp_ref, delta_Cp_ref_val):
    delta_T = temperature - temp_ref
    return delta_Cp_ref_val * delta_T

def delta_TS(temperature, temp_ref, S_ref, delta_Cp_ref_val):
    T_ratio = temperature / temp_ref
    delta_T = temperature - temp_ref
    integral = delta_Cp_ref_val * math.log(T_ratio)
    
    result = delta_T * S_ref + temperature * integral
    return result

def delta_G(delta_H_val, delta_TS_val, g_ref):
    return g_ref + delta_H_val - delta_TS_val

################### Ksp getter ##########################

def get_ksp(delta_g, temperature):
    R = 8.3145
    exponent = - delta_g / (R * temperature)
    return math.exp(exponent)

################ initial guess ###########################

def initial_guess_getter(Ksp_val, compounds_num):
    return math.pow(Ksp_val, 1.0/compounds_num)

################ activity obj func #######################

def optimize_with_system(System, Ksp_val, val_init, pressure, temperature):
    z_ref = [1, 0, 0]
    def liq_act_from_molal(molal):
        prediction = System.LiquidActivityCoefficientAsymMolal(temperature,pressure,[55.5093,molal,molal],z_ref)
        coeff_arr = prediction.values()
        product_of_activity = 1
        for i in range(len(coeff_arr)):
            if i==0: continue
            product_of_activity *= (coeff_arr[i] * molal)
        obj_func = Ksp_val - product_of_activity
        return obj_func  
    return optimize.newton(liq_act_from_molal, val_init)

################ T_range getter ###########################

def range_with_steps(bound_low, bound_up, step) :
    temp_list = []
    temp_item = bound_low

    while temp_item < bound_up :
        temp_list.append(temp_item)
        temp_item = temp_item + step
    
    temp_list.append(bound_up)

    return temp_list

############### csv recorder ##############################

def csv_recorder(file, temperature, solubility):
    file.write(str(temperature) + ',' + str(solubility) + '\n')

############### actual dictionary #########################

data_list = dict({
    'NaCl': {
        'gf_ref': -384.138,
        'Cpo_ref': 0.0505,
        'hf_ref': -411.153,
    },
    'Na+': {
        'gf_ref': -261.905,
        'Cpo_ref': 0.0464,
        'hf_ref': -240.12,
    },
    'Cl-': {
        'gf_ref': -131.228,
        'Cpo_ref': -0.1364,
        'hf_ref': -167.159,
    },
    'KCl': {
        'gf_ref': -409.14,
        'Cpo_ref': 0.0513,
        'hf_ref': -436.744,
    },
    'K+': {
        'gf_ref': -283.27,
        'Cpo_ref': 0.0218,
        'hf_ref': -252.38,
    },
    'K2SO4': {
        'gf_ref': -1321.37,
        'Cpo_ref': 0.13146,
        'hf_ref': -1437.79,
    },
    'SO42-': {
        'gf_ref': -744.53,
        'Cpo_ref': -0.293,
        'hf_ref': -909.27,
    }

})


############### get data from csv #########################

def get_properties(compounds_plus, compounds_minus):
    gf_ref = 0
    Cpo_ref = 0
    hf_ref = 0
    for compound_plus in compounds_plus:
        temp_properties = data_list.get(compound_plus)
        gf_ref += temp_properties.get('gf_ref')
        Cpo_ref += temp_properties.get('Cpo_ref')
        hf_ref += temp_properties.get('hf_ref')
    for compound_minus in compounds_minus:
        temp_properties = data_list.get(compound_minus)
        gf_ref -= temp_properties.get('gf_ref')
        Cpo_ref -= temp_properties.get('Cpo_ref')
        hf_ref -= temp_properties.get('hf_ref')
    return gf_ref * 1000, Cpo_ref * 1000, hf_ref * 1000

