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

phase = "liquid"

#########################################
compounds_toString = ','.join(compounds)
system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)

file_name_temp = "SinglePhaseFugacityCoefficient_" + str(T) + ',' + str(P) + '_' + compounds_toString + "_liquid.csv"

f_csv = open('./predicted_properties/liquid_fugacity_coefficient/' + file_name_temp,"w")

#########################################
x_range = []
m_temp = 0
x_temp = 0
while m_temp <= m_up :
    x_temp = m_temp / ((1000/18.02) + m_temp)
    x_range.append(x_temp)
    m_temp += m_step

print(x_range)
f_csv_firstline = 'x,' + compounds_toString

f_csv.write(f_csv_firstline + '\n')

for x in x_range :
    predicted_property = S.SinglePhaseFugacityCoefficient(T,P,[1-(2*x), x, x],phaseType='liquid')
    print(predicted_property)
    temp_values = predicted_property.values()
    print(temp_values)

    string_to_print = str(x) + ',' + ','.join(map(str, temp_values))

    print(string_to_print)
    f_csv.write(string_to_print + '\n')

print ('\nSUCCESSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS\n\n')

f_csv.close()

