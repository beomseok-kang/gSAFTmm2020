import func_lib as f

compounds = [
    "water",
    "Na+",
    "Cl-"
]

compounds_plus = [
    'Na+',
    'Cl-'
]

compounds_minus = [
    'NaCl'
]

T_low = 280.0
T_up = 500.0
T_step = 1.0

T_ref = 298.15

P = 101325.0

filename = 'Solubility_NaCl.csv'

################ Calculations #############################

delta_g_ref, delta_Cp_ref, delta_H_ref = f.get_properties(compounds_plus, compounds_minus)
S = f.saft_system(compounds)
delta_S_ref = f.S_ref_getter(delta_H_ref, delta_g_ref, T_ref)
T_range = f.range_with_steps(T_low, T_up, T_step)

f_csv = open('./predicted_properties/solubility/' + filename, "w")

for T in T_range:
    delta_g = f.delta_G(
        f.delta_H(T, T_ref, delta_Cp_ref),
        f.delta_TS(T,T_ref,delta_S_ref, delta_Cp_ref),
        delta_g_ref
    )
    Ksp = f.get_ksp(delta_g, T)
    m_initial = f.initial_guess_getter(Ksp, (len(compounds) - 1))
    result = f.optimize_with_system(S, Ksp, m_initial,P, T)
    print((T, result))
    f.csv_recorder(f_csv, T, result)
