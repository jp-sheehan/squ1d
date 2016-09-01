#!/usr/bin/python3
from math import *
import parser
import numpy as np

##### DEFINE HELPER FUNCTIONS #####

# Get next non-empty line of file
def lines(filename):
    with open(filename,'r') as f:
        for l in f:
            c = l.strip()
            if (c):
                if (c == 'EOF'):
                    break;
                else:
                    yield c
    f.close()
    print("done with file " + filename)


def parsedline(iterator,*args):
    parsed = next(iterator).split()
    numparam = len(args)
    if (numparam == 0):
        return
    elif (numparam == 1):
        func = args[0]
        val = func(parsed[0])
        return val
    else:
        vals = []
        for i in range(0,numparam):
            func = args[i]
            vals.append(func(parsed[i]))
        return vals

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def passed(message):
    greencheck = bcolors.OKGREEN + u'\u2713' + bcolors.ENDC
    print(greencheck, message)
    passed.count += 1
    return passed.count

passed.count = 0

def warn(message):
    orangecirc = bcolors.WARNING + '!' + bcolors.ENDC
    print(orangecirc, message)
    warn.count += 1
    return warn.count

warn.count = 0

def err(message):
    redx = bcolors.FAIL + 'X' + bcolors.ENDC
    print(redx, message)
    err.count += 1
    return err.count

err.count = 0

def printreport():
#    print(str(passed.count),"passes,", str(warn.count),\
#            "warnings,", str(err.count),"errors")
    print("--------------------")

    if warn.count == 1:
        ws = ""
    else:
        ws = "s"

    if err.count == 1:
        es = ""
    else:
        es = "s"

    print(str(err.count), "error" + es + ',', \
            str(warn.count), "warning" + ws)


def compileEqn(eqn_str):
    eqn_str = eqn_str.replace('^','**')
    eqn_str = "".join(eqn_str.split()) # remove all whitespace
    if eqn_str[:2] == "if":
        cond, eqn_true, eqn_false = eqn_str[3:-1].split(',')
        if eval(parser.expr(cond).compile()):
            return parser.expr(eqn_true).compile()
        else:
            return parser.expr(eqn_false).compile()
    else:
        return parser.expr(eqn_str).compile()



##### IMPORT FILES #####

# SolverType.inp -> variables
st_iter = lines('SolverType.inp')

solver_type = " ".join(parsedline(st_iter,str,str,str)[1:3])
varname, mesh_read = parsedline(st_iter,str,int)
assert varname == 'MESH_READ', varname + ' should be MESH_READ'
varname, mesh_write = parsedline(st_iter,str,str)
assert varname == 'MESH_WRITE', varname + ' should be MESH_WRITE'
varname, restart = parsedline(st_iter,str,int)
assert varname == 'RESTART', varname + ' should be RESTART'
varname, space_dims, space_symmetry = parsedline(st_iter,str,int,str)
assert varname == 'COORDINATE_TYPE', varname + ' should be COORDINATE_TYPE'
varname, vel_dims = parsedline(st_iter,str,int)
assert varname == 'VECTOR_SIZE', varname + ' should be VECTOR_SIZE'
varname, mesh_length = parsedline(st_iter,str,float)
assert varname == 'MESH_LENGTH', varname + ' should be MESH_LENGTH'
varname, mesh_start = parsedline(st_iter,str,float)
assert varname == 'MESH_START', varname + ' should be MESH_START'
varname, mesh_dims = parsedline(st_iter,str,int)
assert varname == 'MESH_DIMS', varname + ' should be MESH_DIMS'
varname, mesh_weight = parsedline(st_iter,str,float)
assert varname == 'MESH_WEIGHT', varname + ' should be MESH_WEIGHT'
varname, area_over_length = parsedline(st_iter,str,float)
assert varname == 'AREA/LENGTH', varname + ' should be AREA/LENGTH'
varname, ghost_cells = parsedline(st_iter,str,int)
assert varname == 'GHOST_CELLS', varname + ' should be GHOST_CELLS'
varname, Bfield_location, Efield_or_Phi_location = parsedline(st_iter,str,int,int)
assert varname == 'FIELD_LOCATION', varname + ' should be FIELD_LOCATION'
varname, particle_dist = parsedline(st_iter,str,int)
assert varname == 'PARTICLE_DIST', varname + ' should be PARTICLE_DIST'
varname, smoothing = parsedline(st_iter,str,int)
assert varname == 'SMOOTHING', varname + ' should be SMOOTHING'

varname, output_avgs = parsedline(st_iter,str,int)
assert varname == 'OUTPUT_AVGS', varname + ' should be OUTPUT_AVGS'
varname, output_continuum = parsedline(st_iter,str,int)
assert varname == 'OUTPUT_CONTINUUM', varname + ' should be OUTPUT_CONTINUUM'
varname, output_particle = parsedline(st_iter,str,int)
assert varname == 'OUTPUT_PARTICLE', varname + ' should be OUTPUT_PARTICLE'
varname, output_vdf, output_vdf_region = parsedline(st_iter,str,int,int)
assert varname == 'OUTPUT_VDF', varname + ' should be OUTPUT_VDF'
varname, output_rst = parsedline(st_iter,str,int)
assert varname == 'OUTPUT_RST', varname + ' should be OUTPUT_RST'
varname, output_final_numsteps, output_final_skipping = parsedline(st_iter,str,int,int)
assert varname == 'OUTPUT_FINAL', varname + ' should be OUTPUT_FINAL'

varname, length_scale = parsedline(st_iter,str,float)
assert varname == 'LENGTH_SCALE', varname + ' should be LENGTH_SCALE'
varname, time_scale = parsedline(st_iter,str,float)
assert varname == 'TIME_SCALE', varname + ' should be TIME_SCALE'
varname, cfl = parsedline(st_iter,str,float)
assert varname == 'CFL', varname + ' should be CFL'
varname, time_step = parsedline(st_iter,str,float)
assert varname == 'TIME_STEP', varname + ' should be TIME_STEP'
varname, iterations = parsedline(st_iter,str,int)
assert varname == 'ITERATIONS', varname + ' should be ITERATIONS'
varname, print_iterations = parsedline(st_iter,str,int)
assert varname == 'PRINT_ITERATION', varname + ' should be PRINT_ITERATION'
varname, particle_weight = parsedline(st_iter,str, float)
assert varname == 'PARTICLE_WEIGHT', varname + ' should be PARTICLE_WEIGHT'
varname, interpolation_scheme = parsedline(st_iter,str,str)
assert varname == 'INTERPOLATION_SCHEME', varname + ' should be INTERPOLATION_SCHEME'
varname, particle_mover = parsedline(st_iter,str,str)
assert varname == 'PARTICLE_MOVER', varname + ' should be PARTICLE_MOVER'
varname, background_charge = parsedline(st_iter,str,float)
assert varname == 'BACKGROUND_CHARGE', varname + ' should be BACKGROUND_CHARGE'

# SolverInput.inp -> variables
si_iter = lines('SolverInput.inp')

varname, num_species = parsedline(si_iter,str,int)
assert varname == 'SPECIES', varname + ' should be SPECIES'
parsedline(si_iter) #header
species = []
for i in range(0,num_species):
    specie={}
    specie['name'], specie['charge'], specie['mass'], specie['magnetization'] = \
            parsedline(si_iter,str,float,float,int)
    species.append(specie)
assert species[0]['name'] == 'electron', "Electrons should be the first species"

varname, coulomb_collisions = parsedline(si_iter,str,int)
assert varname == 'COULOMB_COLLISIONS', varname + ' should be COULOMB_COLLISIONS'

varname, num_collision_types = parsedline(si_iter,str,int)
assert varname == 'COLLISION_TYPES', varname + ' should be COLLISION_TYPES'
parsedline(si_iter) #header
collisions=[]
for i in range(0,num_collision_types):
    collision={}
    collision['name'], collision['type'], collision['crossSection'] = \
            parsedline(si_iter,str,str,str)
    collisions.append(collision)
varname, n_n = parsedline(si_iter,str,float)
assert varname == 'NeutralDensity', varname + ' should be NeutralDensity'
varname, T_n = parsedline(si_iter,str,float)
assert varname == 'NeutralTemperature', varname + ' should be NeutralTemperature'

varname, num_boundaries = parsedline(si_iter,str,int)
assert varname == 'BOUNDARIES', varname + ' should be BOUNDARIES'
boundaries = []
for i in range(num_boundaries):

    boundary={}
    varname, boundary['name'] = parsedline(si_iter,str,str)
    assert varname == 'Name', varname + ' should be Name'
    varname, boundary['range'] = parsedline(si_iter,str,float)
    assert varname == 'Range', varname + ' should be Range'

    parsedline(si_iter) # header

    bc = {};
    bcname, bc['condition'], bc['value'], bc['distribution'] = \
            parsedline(si_iter,str,str,float,str)
    assert bcname == 'eDensity', bcname + ' should be eDensity'
    boundary[bcname] = bc

    bc = {};
    bcname, bc['condition'], bc['value'], bc['distribution'] = \
            parsedline(si_iter,str,str,float,str)
    assert bcname == 'eTemp', bcname + ' should be eTemp'
    boundary[bcname] = bc

    bc = {};
    bcname, bc['condition'], vx, vy, vz = \
            parsedline(si_iter,str,str,float,float,float)
    assert bcname == 'eVel', bcname + ' should be eVel'
    bc['value'] = [vx, vy, vz]
    boundary[bcname] = bc

    bc = {};
    bcname, bc['condition'], bc['value'], bc['distribution'] = \
            parsedline(si_iter,str,str,float,str)
    assert bcname == 'iDensity', bcname + ' should be iDensity'
    boundary[bcname] = bc

    bc = {};
    bcname, bc['condition'], bc['value'], bc['distribution'] = \
            parsedline(si_iter,str,str,float,str)
    assert bcname == 'iTemp', bcname + ' should be iTemp'
    boundary[bcname] = bc

    bc = {};
    bcname, bc['condition'], vx, vy, vz = \
            parsedline(si_iter,str,str,float,float,float)
    assert bcname == 'iVel', bcname + ' should be iVel'
    bc['value'] = [vx, vy, vz]
    boundary[bcname] = bc

    bc = {};
    bcname, bc['condition'], bc['eqn'] = \
            parsedline(si_iter,str,str,str)
    assert bcname == 'Phi', bcname + ' should be Phi'
    boundary[bcname] = bc

    bc = {};
    bcname, bc['condition'], bx, by, bz = \
            parsedline(si_iter,str,str,str,str,str)
    assert bcname == 'BField', bcname + ' should be BField'
    bc['eqns'] = [bx, by, bz]
    boundary[bcname] = bc

    boundaries.append(boundary)


varname = parsedline(si_iter,str)
assert varname == 'INITIAL_CONDITIONS', varname + " should be INITIAL_CONDITIONS"
initial = {}

ic = {}
icname, ic['value'], ic['distribution'], ic['perturbation'] = \
        parsedline(si_iter,str,float,str,float)
assert icname == 'ElecDensity', icname + ' should be ElecDensity'
initial['eDensity'] = ic

ic = {}
icname, ic['value'], ic['distribution'] = \
        parsedline(si_iter,str,float,str)
assert icname == 'ElecTemperature', icname + ' should be ElecTemperature'
initial['eTemp'] = ic

ic = {}
icname, vx, vy, vz = \
        parsedline(si_iter,str,float,float,float)
assert icname == 'ElecVelocity', icname + ' should be ElecVelocity'
ic['value'] = [vx, vy, vz]
initial['eVel'] = ic

ic = {}
icname, ic['value'], ic['distribution'], ic['perturbation'] = \
        parsedline(si_iter,str,float,str,float)
assert icname == 'IonDensity', icname + ' should be IonDensity'
initial['iDensity'] = ic

ic = {}
icname, ic['value'], ic['distribution'] = \
        parsedline(si_iter,str,float,str)
assert icname == 'IonTemperature', icname + ' should be IonTemperature'
initial['iTemp'] = ic

ic = {}
icname, vx, vy, vz = \
        parsedline(si_iter,str,float,float,float)
assert icname == 'IonVelocity', icname + ' should be IonVelocity'
ic['value'] = [vx, vy, vz]
initial['iVel'] = ic

icname, initial['Phi'] = parsedline(si_iter,str,str)
assert icname == 'Phi', icname + ' should be Phi'

icname, bx, by, bz = parsedline(si_iter,str,str,str,str)
assert icname == 'MagneticField', icname + ' should be MagneticField'
initial['BField'] = [bx, by, bz]

# SPECIAL REGIONS
#varname, numSpecialRegions = parsedline(si_iter, str, int)
#specialRegions = []
#for i in range(numSpecialRegions):
#    region = {}
#
#    varname, region['type'] = parsedline(si_iter, str, str)
#    


##### CALCULATE PARAMETERS #####

# constants

epsilon_0 = 8.854e-12
mu_0 = 1.257e-6
k = 1.381e-23
e = 1.602e-19

# B field interpretation

L_domain = mesh_length - mesh_start
L_cell = L_domain / (mesh_dims - 1)
x_arr = np.linspace(mesh_start,L_domain,mesh_dims)
B_x = np.empty_like(x_arr)
for i in range(x_arr.size):
    x = x_arr[i]
    code_B_x = compileEqn(initial['BField'][0])
    B_x[i] = eval(code_B_x)
B_mag = np.asarray([B_x.min(), B_x.max()])
dB_x = np.abs(np.diff(B_x))
L_B = B_mag.max() / dB_x.max()

A_norm = B_x / B_mag.max()
A_min = B_mag.min()/B_mag.max()

# basic parameters

T_e = initial['eTemp']['value']
T_i = initial['iTemp']['value']
n_e_max = initial['eDensity']['value']
n_e = (A_min * n_e_max, n_e_max)
n_i_max = initial['iDensity']['value']
n_i = (A_min * n_i_max, n_i_max)
m_e = species[0]['mass']
m_i = species[1]['mass']
q_e = species[0]['charge']
q_i = species[1]['charge']
v_e_th = sqrt(8 * k * T_e / (pi * m_e)) # average thermal velocity
v_i_th = sqrt(8 * k * T_i / (pi * m_i))
c_s = sqrt(k * T_e / m_i)
v_i_b = sqrt(log(sqrt(m_i/(2*pi*m_e)))*k*T_e/m_i) # order of magnitude beam velocity

# neutral collisions
# nn and Tn defined above

sigma_in = 1e-18 # general approximation, charge exchange
sigma_en = 1e-20 # general approximation, elastic
mfp_in = 1/(n_n * sigma_in)
mfp_en = 1/(n_n * sigma_en)
nu_in = v_i_b/mfp_in
nu_en = v_e_th/mfp_en
tau_in = 1/nu_in
tau_en = 1/nu_en

# simulation size

vol = 1 * 1 * mesh_length
N_tot_i = initial['iDensity']['value']/particle_weight * vol
N_tot_e = initial['eDensity']['value']/particle_weight * vol

N_per_cell_e = N_tot_e/mesh_dims
N_per_cell_i = N_tot_i/mesh_dims

T_tot = 0.01 * iterations * N_tot_e

# plasma properties

lambda_de = np.sqrt(np.divide(epsilon_0 * k * T_e, np.multiply(n_e,q_e**2)))
lambda_di = np.sqrt(np.divide(epsilon_0 * k * T_i, np.multiply(n_i,q_i**2)))

omega_pe = np.sqrt(np.divide(np.multiply(n_e, q_e**2), m_e * epsilon_0))
omega_pi = np.sqrt(np.divide(np.multiply(n_i, q_i**2), m_i * epsilon_0))

Lambda = np.multiply(4/3 * pi, np.multiply(n_e, np.power(lambda_de,3)))

omega_ce = np.abs(np.multiply(B_mag,q_e/m_e))
omega_ci = np.abs(np.multiply(B_mag,q_i/m_i))

r_Le = np.divide(v_e_th,omega_ce)
r_Li = np.divide(v_e_th,omega_ci)

r_L_inertial = np.add(r_Le.max(), np.multiply(r_Li.max(), v_i_th/v_e_th))
r_L_hybrid = np.sqrt(np.multiply(r_Le.max(), r_Li.max()))

# time scales

tau_i_cell = L_cell / v_i_b
tau_e_cell = L_cell / v_e_th
tau_i_domain = L_domain / v_i_b
tau_e_domain = L_domain / v_e_th

# file sizes

# bytes per line, approximate
bpl_cpField = 230
bpl_particles = 50
bpl_phasespace = 40
bpl_vdf = 75
bpl_outField = 100
bpl_energy = 75

# lines per file, approximate
lpf_cpField = mesh_dims
lpf_particles = N_tot_e
lpf_phasespace = 100 * mesh_dims
lpf_vdf = 100
lpf_outField = mesh_dims
lpf_energy = iterations

if output_vdf_region == -1:
    vdf_frac = 1
else:
    vdf_frac = 1/mesh_dims

sizePerWrite = num_species * (bpl_cpField * lpf_cpField + \
        output_particle * bpl_particles * lpf_particles + \
        output_continuum * bpl_phasespace * lpf_phasespace + \
        output_vdf * vdf_frac * bpl_vdf * lpf_vdf) + \
        output_continuum * bpl_outField * lpf_outField
N_writes = iterations//print_iterations + \
        output_final_numsteps//output_final_skipping
writeSizeTot = sizePerWrite * N_writes + output_avgs * bpl_energy * lpf_energy



##### CHECK PARAMETERS #####

print("---Checking unused parmeters---")

if mesh_weight != 1: warn("MESH_WEIGHT will be unused")
if length_scale != 1: warn("LENGTH_SCALE will be unused")
if time_scale != 1: warn("TIME_SCALE is will be unused")
if cfl != 1: warn("CFL will be unused")

N_we = warn.count + err.count
if N_we == 0: passed("Passed")

print("---Checking for consistency---")
if not(boundaries[0]['BField']['eqns'] == boundaries[1]['BField']['eqns']  == initial['BField']):
    err("Inconsistent magnetic field definitions")
if not(boundaries[0]['Phi']['eqn']  == boundaries[1]['Phi']['eqn']  == initial['Phi'] ):
    err("Inconsistent potential definitions")
if not(boundaries[0]['eDensity']['value'] == boundaries[1]['eDensity']['value'] == initial['eDensity']['value']):
    err("Inconsistent electron density definitions")
if not(boundaries[0]['iDensity']['value'] == boundaries[1]['iDensity']['value'] == initial['iDensity']['value']):
    err("Inconsistent ion density definitions")
if initial['iDensity']['value'] != initial['eDensity']['value']: err("Initial conditions are non-quasineutral")

if warn.count + err.count == N_we: passed("Passed")
N_we = warn.count + err.count

print("---Checking for uncommon settings---")

if mesh_write != "TECPLOT": err("Use TECPLOT file type for postprocessing compatibility")
if vel_dims != 3: err("Use 3 velocity dimensions")
if Bfield_location != 1: err("Bfield should be defined on corners")
if Efield_or_Phi_location != 1: err("Efield or Phi should be defined on corners")
if not(initial['BField'][1] == initial['BField'][2] == '0.0'): warn('Simulation space does not follow mangetic field line.  Checker may miss some things.')

if warn.count + err.count == N_we: passed("Passed")
N_we = warn.count + err.count

print("---Checking simulation parameters---")

omega_str = u"\u03c9"
Delta_str = u"\u0394"
lambda_str = u"\u03bb"
nu_str = u"\u039d"

# time step
tFactor_e = omega_pe.max() * time_step
if tFactor_e > 0.2: err("Time step too long  (" + omega_str + "_pe*" + Delta_str + "t=%.2e)" % tFactor_e)
tFactor_i = omega_pi.max() * time_step
if tFactor_i > 0.2: err("Time step too long  (" + omega_str + "_pi*" + Delta_str + "t=%.2e)" % tFactor_i)

# simulation time
TFactor_e = (time_step * iterations) / tau_e_domain # transit time / simulation time
if TFactor_e < 1: err("Not enough time for electrons to transit entire domain (tau_e/N" + Delta_str + "t=%.2e)" % TFactor_e)
TFactor_i = (time_step * iterations) / tau_i_domain
if TFactor_i < 1: warn("Not enough time for ions to transit entire domain (tau_i/N" + Delta_str + "t=%.2e)" % TFactor_i)

# cell size
xFactor_e = L_cell / lambda_de.min()
if xFactor_e > 0.5: err("Cell size too large (" + Delta_str + "x/" + lambda_str + "_de=%.2e)" % xFactor_e)
xFactor_i = L_cell / lambda_di.min()
if xFactor_i > 0.5: err("Cell size too large (" + Delta_str + "x/" + lambda_str + "_di=%.2e)" % xFactor_i)

# particles per cell
ppc_e = N_tot_e // (mesh_dims-1)
ppc_i = N_tot_i // (mesh_dims-1)
if (ppc_e < 400) | (ppc_i < 400): err("Too few particles per cell (%d electrons, %d ions)" % (ppc_e, ppc_i))

T_days = T_tot / (1000*60*60*24)
if T_days > 1: warn("Simulation will take ~%d days on a single processor" % T_days)

if warn.count + err.count == N_we: passed("Passed")
N_we = warn.count + err.count

print("---Checking plasma parameters---")

wFactor_e = nu_en/omega_ce.min()
if wFactor_e > 1: warn("Neutral collisions will prevent electron confinement (" + nu_str + "_en/" + omega_str + "_ce=%.2e" % wFactor_e)

riFactor = r_L_inertial.max() / L_B
if riFactor > 1e-2: err("Inertial Larmor radius is too large (r_L/L_B = %.2e)" % riFactor)

if warn.count + err.count == N_we: passed("Passed")
N_we = warn.count + err.count

print("---Checking file sizes---")

if writeSizeTot > 1e12: err("Over 1TB will be written")
elif writeSizeTot > 1e11: err("Over 100GB will be written")
elif writeSizeTot > 1e10: warn("Over 10GB will be written")
elif writeSizeTot > 1e9: warn("Over 1GB will be written")
elif writeSizeTot > 1e8: warn("Over 100MB will be written")

if warn.count + err.count == N_we: passed("Passed")
N_we = warn.count + err.count



printreport()
