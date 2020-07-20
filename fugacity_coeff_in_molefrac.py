import gSAFTmm

compounds = [
    "water",
    "ethanol_2"
]

T = 298
P = 101325

m_low = 0.0
m_up = 1.0
m_step = 0.01

phase = "liquid"

#########################################
compounds_toString = ','.join(compounds)
system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)

file_name_temp = "SinglePhaseFugacityCoefficient_" + str(T) + ',' + str(P) + '_' + compounds_toString + "_liquid.csv"

f_csv = open('./predicted_properties/liquid_fugacity_coefficient/' + file_name_temp,"w")

#########################################
m_range = []
m_temp = m_low
while m_temp <= m_up :
    m_range.append(m_temp)
    m_temp += m_step
m_range.append(m_up)

f_csv_firstline = 'm,' + compounds_toString

f_csv.write(f_csv_firstline + '\n')

for m in m_range :
    predicted_property = S.SinglePhaseFugacityCoefficient(T,P,[1-m, m],phaseType='liquid')
    temp_values = predicted_property.values()

    string_to_print = str(m) + ',' + ','.join(map(str, temp_values))

    print(string_to_print)
    f_csv.write(string_to_print + '\n')

print ('\nSUCCESSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS\n\n')

f_csv.close()

