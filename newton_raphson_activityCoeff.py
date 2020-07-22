import func_lib as f

compounds = [
    "water",
    "Na+",
    "Cl-"
]

T_low = 280.0
T_up = 500.0
T_step = 1.0

T_ref = 298.15

P = 101325.0

filename = 'Solubility_NaCl.csv'

Nagf_ref = -261.905
Clgf_ref = -131.228
NaClgf_ref = -384.138

Na_Cpo_ref = 0.0464
Cl_Cpo_ref = -0.1364
NaCl_Cpo_ref = 0.05050

Nahf_ref = -240.120
Clhf_ref = -167.159
NaClhf_ref = -411.153

delta_g_ref = (Nagf_ref + Clgf_ref - NaClgf_ref) * 1000.0
delta_Cp_ref = (Na_Cpo_ref + Cl_Cpo_ref - NaCl_Cpo_ref) * 1000.0
delta_H_ref = (Nahf_ref + Clhf_ref - NaClhf_ref) * 1000.0

################ Calculations #############################

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
