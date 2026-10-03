"""one Recola session: V(Delta=0)/B for a given m_b, nf scheme, scale, momenta"""
import sys, math, json
sys.path.append("/home/cbiello/Packages/recola-heft-py/lib/python3.13/site-packages")
from pyrecola import *
cfg = json.loads(sys.argv[1])
mb, nfs, mu, aS, proc = cfg["mb"], cfg["nf"], cfg["mu"], cfg["aS"], cfg["proc"]
p = cfg["p"]
set_output_file_rcl('amp.rcl')
set_print_level_squared_amplitude_rcl(0); set_print_level_parameters_rcl(0)
set_mu_uv_rcl(mu); set_mu_ir_rcl(mu); set_mu_ms_rcl(mu)
set_delta_ir_rcl(0.,0.); set_dynamic_settings_rcl(1)
set_parameter_rcl('WW',0.); set_parameter_rcl('WZ',0.)
set_parameter_rcl('WH',0.); set_parameter_rcl('WT',0.); set_parameter_rcl('WB',0.)
set_parameter_rcl('MB',cfg.get('MBloop',mb)); set_parameter_rcl('ymb',0.)
set_parameter_rcl('MT',173.2); set_parameter_rcl("MH",cfg.get("MH",125.))
set_parameter_rcl('cgghfin',0.); set_parameter_rcl('aS',aS)
set_renoscheme_rcl('dZgs_QCD2','Nf%d'%nfs)
define_process_rcl(1,proc,'NLO'); generate_processes_rcl()
out={}
for lbl,(d1,d2) in {"f":(0.,0.),"p1":(1.,0.),"p2":(0.,1.)}.items():
    set_mu_uv_rcl(mu); set_mu_ir_rcl(mu); set_mu_ms_rcl(mu); set_delta_ir_rcl(d1,d2)
    compute_process_rcl(1,p,'NLO')
    B = get_squared_amplitude_rcl(1,'LO', pow=[8,2])
    V = get_squared_amplitude_rcl(1,'NLO',pow=[10,2])
    out[lbl]=V/B/(aS/2/math.pi); out["B"]=B
print(json.dumps(out))
