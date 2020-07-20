import gSAFTmm

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

Ksp = 39.51

type_activity_coefficient = "asymmetric" # either symmetric or asymmetric

#################################################################################

z_ref = []
if type_activity_coefficient == "symmetric" :
    z_ref = [0.001, 1, 1]
else :
    z_ref = [55.51, 0, 0]

#################################################################################

compounds_toString = ','.join(compounds)
system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)

file_name_temp = "activityCoeff_" + type_activity_coefficient + str(T) + ',' + str(P) + '_' + compounds_toString + ".csv"

f_csv = open('./predicted_properties/liquid_activity_coefficient_fromgSAFT/' + file_name_temp,"w")

#################################################################################

m_range = []
m_temp = m_low
while m_temp <= m_up :
    m_range.append(m_temp)
    m_temp += m_step

f_csv_firstline = 'm,'
for compound in compounds :
    f_csv_firstline += compound + ','
f_csv_firstline += type_activity_coefficient

f_csv.write(f_csv_firstline + '\n')

#################################################################################

fugacity_coeff_ref = S.SinglePhaseFugacityCoefficient(T,P,z_ref,phaseType='liquid')
fugacity_coeff_ref_values = fugacity_coeff_ref.values()
print(fugacity_coeff_ref_values)

for m in m_range :
    predicted_property = S.SinglePhaseFugacityCoefficient(T,P,[55.51,m,m],phaseType='liquid')
    # the predicted property will be given in dictionary form
    temp_fugacity_coeff_values = predicted_property.values()
    print(temp_fugacity_coeff_values)
    temp_activity_coeff_values = []
    for i in range(len(temp_fugacity_coeff_values)) :
        temp_activity_coeff = temp_fugacity_coeff_values[i] / fugacity_coeff_ref_values[i]
        temp_activity_coeff_values.append(temp_activity_coeff)
    print(temp_activity_coeff_values)
    string_to_print = str(m) + ',' + ','.join(map(str, temp_activity_coeff_values))
    print(string_to_print)
    f_csv.write(string_to_print + '\n')

print('\nsuccess!\n\n')

f_csv.close()