#!/usr/bin/env python
"""
g g -> b~ b H, with HEFT_2.2.3 ("SM (QCD) + HEFT").

Remember:

  export LD_LIBRARY_PATH=/home/colomba/programs/recola2-collier-2.2.4/COLLIER-1.2.5:/home/colomba/programs/recola2-collier-2.2.4/model_files_2.2.4/HEFT_2.2.3:/home/colomba/programs/recola2-collier-2.2.4-heft-install/lib

  /usr/bin/python3 recola_Hbb.py

"""

from __future__ import print_function
import sys
import math

if not sys.executable.startswith("/usr/bin/python3"):
    sys.exit(
        "Questo script va lanciato con /usr/bin/python3, non con "
        + sys.executable
        + " -- vedi il commento in testa al file."
    )

sys.path.append("/home/colomba/programs/recola2-collier-2.2.4-heft-install/lib/python3.8/site-packages")
# local shared Recola2 build with the python interface (yueh)
sys.path.append("/home/cbiello/Packages/recola-heft-py/lib/python3.13/site-packages")

from pyrecola import *

set_output_file_rcl('*')
set_print_level_squared_amplitude_rcl(2)
set_print_level_parameters_rcl(2)

set_mu_uv_rcl(173.2)
set_mu_ir_rcl(173.2)
set_mu_ms_rcl(173.2)
#use_dim_reg_soft_rcl()
set_delta_ir_rcl(0., math.pi**2/6.)
# lets delta_IR / mu be changed again after generate_processes_rcl,
# needed for the finite-remainder conversion to Bayu
set_dynamic_settings_rcl(1)
#set_on_shell_scheme_rcl()

set_parameter_rcl('WW', 0.)
set_parameter_rcl('WZ', 0.)
#set_parameter_rcl('MB', 5.)
set_parameter_rcl('ymb', 0.)
set_parameter_rcl('yb', 0.)
#set_parameter_rcl("cggh", 2.)

# ---------------------------------------------------------------------
# Process
# ---------------------------------------------------------------------
proc_id = 1

define_process_rcl(proc_id, 'g g -> b~ b H', 'NLO')

generate_processes_rcl()

# with dynamic settings on, the scales must be re-applied after generation
set_mu_uv_rcl(173.2)
set_mu_ir_rcl(173.2)
set_mu_ms_rcl(173.2)
set_delta_ir_rcl(0., math.pi**2/6.)

# ---------------------------------------------------------------------
# Kinematic point
# ---------------------------------------------------------------------
p1_g    = [500.00000000000000000,      0.00000000000000000,      0.00000000000000000,    500.00000000000000000]
p2_g    = [500.00000000000000000,      0.00000000000000000,      0.00000000000000000,   -500.00000000000000000]
p3_bbar = [439.78021928785585715,    162.50680658487215169,    364.09047948353401125,   -185.57020730686238608]
p4_b    = [349.14239937112841972,    -17.57847102882445967,   -333.45082223855939674,    101.99000707591839898]
p5_H    = [211.07738134101575156,   -144.92833555604772755,    -30.63965724497451859,     83.58020023094401552]

# ---------------------------------------------------------------------

psp = [p1_g, p2_g, p3_bbar, p4_b, p5_H]

compute_process_rcl(proc_id, psp, 'NLO')
#set_qcd_rescaling_rcl(True)
#compute_running_alphas_rcl(173.2,5,2)

amp_LO = get_squared_amplitude_rcl(proc_id, 'LO', pow=[8, 2])
print("LO,  |A0|^2 =", amp_LO)

#get_amplitude_rcl(1, 'LO', [1,2,0,4,0], [-1,-1,+1,-1,0], pow=[8,2])

amp_NLO_0 = get_squared_amplitude_rcl(proc_id, 'NLO', pow=[10, 2])
print("NLO =", amp_NLO_0)
print("K-factor (NLO/LO) =", amp_NLO_0 / amp_LO if amp_LO != 0 else float('nan'))

# ---------------------------------------------------------------------
# Finite remainder in Bayu's convention
#
# Recola returns  V = c2*Delta_IR2 + c1*Delta_IR + f , so re-evaluating at
# Delta_IR2 = 0 and 1 gives the eps^0 coefficient f and the double pole c2.
# Two conventions then separate Recola from Bayu:
#
#   1. eps-expansion: the paper factors out S_eps = (4pi)^eps e^(-eps gammaE),
#      Recola's Delta_IR2 = zeta2 corresponds to (4pi)^eps/Gamma(1-eps).
#      Since e^(eps gE)/Gamma(1-eps) = 1 + zeta2 eps^2/2 + ..., the eps^0
#      parts differ by the usual annoying c2*zeta2/2.
#   2. Wilson coefficient: Bayu's factor C1 out of the amplitude, while
#      the Recola HEFT model has dcgghfin_QCD2 = 11 aS/(4 pi) as we saw
#      in the meeting (it does not affect at mu=mt the ration V/B).
#
# The paper's IR subtraction is SCET so no shift of finite part as in Catani..
# ---------------------------------------------------------------------

aS   = 0.118
unit = aS/(2.*math.pi)  # normalisation of V/B
zeta2 = math.pi**2/6.

def virt(d2):
    set_mu_uv_rcl(173.2); set_mu_ir_rcl(173.2); set_mu_ms_rcl(173.2)
    set_delta_ir_rcl(0., d2) #we compute the virtual with different eps NORMALISATIONS
    compute_process_rcl(proc_id, psp, 'NLO')
    return get_squared_amplitude_rcl(proc_id, 'NLO', pow=[10, 2])

f_eps0 = virt(0.) / amp_LO / unit          # eps^0 coefficient
c2     = virt(1.) / amp_LO / unit - f_eps0 # 1/eps^2 coefficient = -sum_i C_i
finrem = f_eps0 + c2*zeta2/2. - 2.*(11./4.)*(aS/math.pi)/unit # the last remove the wired Recola wilson normalisation

print()
print("V/B, BLHA (Delta_IR2 = zeta2)      =", amp_NLO_0/amp_LO/unit)
print("eps^0 coefficient f                =", f_eps0)
print("1/eps^2 coefficient c2 (-sum C_i)  =", c2)
print("+ c2*zeta2/2   (S_eps convention)  =", f_eps0 + c2*zeta2/2.)
print("- C1 O(aS) piece  -> finite rem.   =", finrem)
print("  same, as Msq1L/Msq0L             =", finrem*unit)
print("  Bayu reference                   =", 0.19644658463793435)

reset_recola_rcl()
