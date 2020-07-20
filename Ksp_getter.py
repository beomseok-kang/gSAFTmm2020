import gSAFTmm
import math

compounds = [
    "water",
    "Na+",
    "Cl-"
]

T = 298
P = 101325

m_low = 0.0
m_up = 1.0
m_step = 0.01

g_form_NaCl = -384024
g_form_Na = 123123
g_form_Cl = 123123

##########################################################################################

compounds_toString = ','.join(compounds)
system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)

file_name_temp = "Ksp_" + str(T) + ',' + str(P) + '_' + compounds_toString + ".csv"

f_csv = open('./predicted_properties/Ksp/' + file_name_temp,"w")

##########################################################################################

## m_range (molality range)

m_range = []
m_temp = m_low
while m_temp <= m_up :
    m_range.append(m_temp)
    m_temp += m_step
m_range.append(m_up)

## first line

f_csv_firstline = 'm,'
for compound in compounds :
    f_csv_firstline += 'chemical potential(' + compound + '),'
f_csv_firstline += 'Ksp'

f_csv.write(f_csv_firstline + '\n')

##########################################################################################

for m in m_range :
    x_i = m / ( 2 * m + 55.51)
    x_water = 1 - 2 * x_i
    prediction = S.SinglePhaseChemicalPotential(T,P,[x_water, x_i, x_i], phaseType="liquid")

    prediction_to_list = prediction.values()

    ###### Ksp Calculation ######
    exponent = - (prediction_to_list[1] + prediction_to_list[2] - g_form_NaCl) / (8.314 * T)
    Ksp = math.exp(exponent)
    #############################

    string_to_print = str(m) + ',' + ','.join(map(str, prediction_to_list)) + ',' + str(Ksp)

    print(string_to_print)
    f_csv.write(string_to_print + '\n')

print('finished')

f_csv.close()

