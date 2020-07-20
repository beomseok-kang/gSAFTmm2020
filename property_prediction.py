import gSAFTmm

#### compounds ######################################################################################

compounds = [
    "water"
    , 
    "Na+"
    ,
    "Cl-"
]

#### Inputs #########################################################################################

# property to predict
property_to_predict = "SinglePhaseFugacityCoefficient"

# T upper and lower limit, and the step in [K]
T_low = 298
T_up = 328
T_step = 1

# T upper and lower limit, and the step in [Pa]
P_low = 101325
P_up = 131325
P_step = 1000

# compounds mole fraction (does not necessarily have to sum up to 1, though)
compounds_mole_frac = [
    # First compound
    0.9
    , # Second compound
    0.05
    , # ...so on
    0.05
]

# phase: either vapour, liquid, or leave blank(most_stable)
phase = "liquid"






#####################################################################################################
#####################################################################################################
#### SAFT System (Don't edit the codes below!) ######################################################
#####################################################################################################
#####################################################################################################

compounds_toString = ','.join(compounds)

system_string = "GC_Mie_Databank_TL_MW_050620.xml<" + compounds_toString + ">"

S = gSAFTmm.System(system_string)


#### What to measure ################################################################################

def range_with_steps(bound_low, bound_up, step) :
    temp_list = []
    temp_item = bound_low

    while temp_item < bound_up :
        temp_list.append(temp_item)
        temp_item = temp_item + step
    
    temp_list.append(bound_up)

    return temp_list
        
#### save as csv #####################################################################################

file_name_temp = property_to_predict + "_" + compounds_toString + "_" + str(compounds_mole_frac) + ("_" + phase if phase else "") + ".csv"

f_csv = open('./predicted_properties/' + file_name_temp,"w")

#### property specification ##########################################################################

property_type = ""
pure_component_vap_liq_equilibrium = [
    "BoilingPoint",
    "VapourPressure",
    "SaturationLiquidDensity",
    "SaturationVapourDensity",
    "VaporizationEnthalpy"
]
single_phase_properties = [
    "SinglePhaseChemicalPotential",
    "SinglePhaseCompressibilityFactor",
    "SinglePhaseCpCv",
    "SinglePhaseDensity",
    "SinglePhaseEnergy",
    "SinglePhaseEnthalpy",
    "SinglePhaseEntropy",
    "SinglePhaseExcessEnthalpy",
    "SinglePhaseExcessGibbsEnergy",
    "SinglePhaseFugacityCoefficient",
    "SinglePhaseFugacityCoefficientLn",
    "SinglePhaseGibbsFreeEnergy",
    "SinglePhaseHeatCapacity",
    "SinglePhasePartialEnthalpy",
    "SinglePhaseSpeedOfSound",
    "SinglePhaseVolume"
]
phase_boundary_methods_with_T = [
    "BubblePressure",
    "BubblePressureCompositions",
    "DewPressure",
    "DewPressureCompositions"
]
phase_boundary_methods_with_P = [
    "BubbleTemperature",
    "BubbleTemperatureCompositions",
    "DewTemperature",
    "DewTemperatureCompositions"
]
total_properties = [
    "CompressibilityFactor",
    "CpCv",
    "Density",
    "Energy",
    "Enthalpy",
    "Entropy",
    "GibbsEnergy",
    "HeatCapacity",
    "SpeedOfSound",
    "Volume"
]
liquid_phase_properties = [
    "LiquidActivityCoefficient",
    "LiquidActivityCoefficientLn",
    "LiquidpH"
]

if property_to_predict in pure_component_vap_liq_equilibrium :
    property_type = "pure_component_vap_liq_equilibrium"
elif property_to_predict in single_phase_properties :
    property_type = "single_phase_properties"
elif property_to_predict in phase_boundary_methods_with_T :
    property_type = "phase_boundary_methods_with_T"
elif property_to_predict in phase_boundary_methods_with_P :
    property_type = "phase_boundary_methods_with_P"
elif property_to_predict in total_properties :
    property_type = "total_properties"
elif property_to_predict in liquid_phase_properties:
    property_type = "liquid_phase_properties"
else :
    property_type = "invalid type or flash properties or Liquid Phase AsymMolal or AsymMolar properties"
    print(property_type)

#### actual calculations #############################################################################
success_message = """

      .--..--..--..--..--..--. 
    .' \  (`._   (_)     _   /      __________________________________
  .'    |  '._)         (_)  |     |                                  |
  \ _.')\      .----..---.   /     |   All Calculations Successful!   |
  |(_.'  |    /    .-\-.  \  |    /___________________________________|
  \     0|    |   ( O| O) | o|
   |  _  |  .--.____.'._.-.  |
   \ (_) | o         -` .-`  |
    |    \   |`-._ _ _ _ _\ /
    \    |   |  `. |_||_|   |
    | o  |    \_      \     |     -.   .-.
    |.-.  \     `--..-'   O |     `.`-' .'
  _.'  .' |     `-.-'      /-.__   ' .-'
.' `-.` '.|='=.='=.='=.='=|._/_ `-'.'
`-._  `.  |________/\_____|    `-.'
   .'   ).| '=' '='\/ '=' |
   `._.`  '---------------'
           //___\   //___|
             ||       ||
             ||_.-.   ||_.-.
            (_.--__) (_.--__)

"""

if T_low > T_up or P_low > P_up :
    print("Error: T_low should not be higher than T_up and P_low should not be higher than P_up")
    
elif property_type == "pure_component_vap_liq_equilibrium" :
    if property_to_predict == "BoilingPoint" :
        f_csv.write(property_to_predict + '\n')
        predicted_property = getattr(S, property_to_predict)()

        string_to_print = str(predicted_property)

        print(string_to_print)
        f_csv.write(string_to_print + '\n')

        print(success_message)

    else :
        f_csv.write('T[K],' + property_to_predict + '\n')

        for T in range_with_steps(T_low, T_up, T_step):
            predicted_property = getattr(S, property_to_predict)(T)
            
            string_to_print = str(T) + "," + str(predicted_property)
            
            print(string_to_print)
            f_csv.write(string_to_print + '\n')

        print(success_message)

elif property_type == "single_phase_properties" :
    phase = phase if phase else "most_stable"
    f_csv.write('T[K],P[kPa]' + property_to_predict + ',phase:' + phase + '\n')

    for T in range_with_steps(T_low, T_up, T_step):
        for P in range_with_steps(P_low, P_up, P_step) :
            predicted_property = getattr(S, property_to_predict)(T, P, compounds_mole_frac, phaseType=phase)
            
            string_to_print = str(T) + "," + str(P) + "," + str(predicted_property)
            
            print(string_to_print)
            f_csv.write(string_to_print + '\n')
    print(success_message)

elif property_type == "phase_boundary_methods_with_T" :
    f_csv.write('T[K],' + property_to_predict + '\n')

    for T in range_with_steps(T_low, T_up, T_step):
        predicted_property = getattr(S, property_to_predict)(T, compounds_mole_frac)
        
        string_to_print = str(T) + "," + str(predicted_property)
        
        print(string_to_print)
        f_csv.write(string_to_print + '\n')
    print(success_message)

elif property_type == "phase_boundary_methods_with_P" :
    f_csv.write('P[kPa],' + property_to_predict + '\n')

    for P in range_with_steps(P_low, P_up, P_step):
        predicted_property = getattr(S, property_to_predict)(P, compounds_mole_frac)
        
        string_to_print = str(P) + "," + str(predicted_property)
        
        print(string_to_print)
        f_csv.write(string_to_print + '\n')
    print(success_message)

elif property_type == "total_properties" :
    f_csv.write('T[K],P[kPa]' + property_to_predict + '\n')

    for T in range_with_steps(T_low, T_up, T_step):
        for P in range_with_steps(P_low, P_up, P_step) :
            predicted_property = getattr(S, property_to_predict)(T, P, compounds_mole_frac)
            
            string_to_print = str(T) + "," + str(P) + "," + str(predicted_property)
            
            print(string_to_print)
            f_csv.write(string_to_print + '\n')
    print(success_message)

elif property_type == "liquid_phase_properties" :
    f_csv.write('T[K],P[kPa]' + property_to_predict + '\n')

    for T in range_with_steps(T_low, T_up, T_step):
        for P in range_with_steps(P_low, P_up, P_step) :
            predicted_property = getattr(S, property_to_predict)(T, P, compounds_mole_frac)
            
            string_to_print = str(T) + "," + str(P) + "," + str(predicted_property)
            
            print(string_to_print)
            f_csv.write(string_to_print + '\n')
    print(success_message)

else :
    print('ERROR - Invalid method: please re-check the property to measure.')

#### all tasks done ###################################################################################

f_csv.close()

