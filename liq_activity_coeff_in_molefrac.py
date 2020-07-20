import gSAFTmm

compounds = [
    "water",
    "ethanol_2"
]

T = 298

P = 101325

# [M] , ethanol
x_low = 0.0
x_up = 1.0
x_step = 0.01

#######################################
z = [55.348, 0.0, 0.0]
# reference molarity of water is 55.348
#######################################

z_ref = [1.0 , 0]
# pure water

############################################

compounds_toString = ','.join(compounds)
system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)

file_name_temp = "LiqActivityCoefficient_" + str(T) + ',' + str(P) + '_' + compounds_toString + "_" + str(z) + ".csv"

f_csv = open('./predicted_properties/liquid_activity_coefficient_molar/' + file_name_temp,"w")

############################################

x_range = []
x_temp = x_low
while x_temp <= x_up :
    x_range.append(x_temp)
    x_temp += x_step
x_range.append(x_up)

f_csv_firstline = 'x,'
for compound in compounds :
    f_csv_firstline += compound + ','

f_csv.write(f_csv_firstline + '\n')

for x in x_range :
    predicted_property = S.LiquidActivityCoefficient(T,P,[1-x,x])
    # the predicted property will be given in dictionary form
    temp_values = predicted_property.values()

    string_to_print = str(x) + ','
    for value in temp_values :
        string_to_print += str(value) + ','

    print(string_to_print)
    f_csv.write(string_to_print + '\n')
print('\nSUCCESSSSSSSSSSSSSSSSSSSSSSSSSSSSSS\n\n')

f_csv.close()