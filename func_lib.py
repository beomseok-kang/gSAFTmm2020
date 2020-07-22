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

################ initial guess ############################

def initial_guess_getter(Ksp_val, compounds_num):
    return math.exp((1.0/compounds_num) * math.log(Ksp_val))

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
