import gSAFTmm
import math

compounds = [
    "water",
    "Na+",
    "Cl-"
]

T = 298

P = 101325
Ksp = 39.51

# [M]
m_low = 0.0
m_up = 7.0
m_step = 0.05

#######################################
z = [55.402, 0.0, 0.0]
# reference molality of water is 55.402
#######################################

z_ref = [55.402 , 0, 0]
# pure water

############################################

compounds_toString = ','.join(compounds)
system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)

file_name_temp = "LiqActivityCoefficientAsymMolal_" + str(T) + ',' + str(P) + '_' + compounds_toString + "_" + str(z) + ".csv"

f_csv = open('./predicted_properties/liquid_activity_coefficient_molal/' + file_name_temp,"w")

############################################

m_range = []
m_temp = m_low
while m_temp <= m_up :
    m_range.append(m_temp)
    m_temp += m_step
m_range.append(m_up)

f_csv_firstline = 'm,'
for compound in compounds :
    f_csv_firstline += compound + ','
f_csv_firstline += 'mole fraction(ion),solubility (molal)'
f_csv.write(f_csv_firstline + '\n')

for m in m_range :
    string_to_print = str(m) + ','

    predicted_property = S.LiquidActivityCoefficientAsymMolal(T,P,[55.402,m,m],z_ref)
    asym_molar = S.LiquidActivityCoefficientAsymMolar(T,P,[55.402,m,m], z_ref)

    # the predicted property will be given in dictionary form
    temp_activity_coeff_values = predicted_property.values()
    temp_asym_molar_activity_coeff = asym_molar.values()

    # solubility_value = Ksp
    for i in range(len(temp_activity_coeff_values)) :
        string_to_print += str(temp_activity_coeff_values[i]) + ','
        string_to_print += str(temp_asym_molar_activity_coeff[i])+','
        if i == 0 :
            continue
    #     solubility_value = solubility_value / temp_activity_coeff_values[i]
    
    mole_frac = m / (2 * m + 55.402)

    string_to_print += str(mole_frac) + ','
    # string_to_print += str(solubility_value)

    print(string_to_print)
    f_csv.write(string_to_print + '\n')
print('\nSUCCESSSSSSSSSSSSSSSSSSSSSSSSSSSSSS\n\n')

f_csv.close()