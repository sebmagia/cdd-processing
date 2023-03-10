import glob, re, os

import numpy as np
import matplotlib.pyplot as plt

import swprepost

def plot_target(target):
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6, 3), dpi=150)
    target.plot(x="frequency", y="velocity", ax=axs[0])
    target.plot(x="wavelength", y="velocity", ax=axs[1])
    axs[1].set_ylabel("")
    axs[1].legend()
    return (fig, axs)

print("Imports successful, you may proceed.")

# Approach 1: Import from comma seperated text file (see swprepost documentation for details).
target = swprepost.Target.from_csv("cdd00.csv")
name="/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO/cdd00.csv"
target = swprepost.Target.from_csv("/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO/cdd00.csv")

# Approach 2: Import from version 2.X.X dinver-style text file (see swprepost documentation for details).
# target = swprepost.Target.from_txt_dinver("example_dv2.txt", version="2")

# Approach 3: Import from version 3.X.X dinver-style text file (see swprepost documentation for details).
# target = swprepost.Target.from_txt_dinver("example_dv3.txt", version="3")


fig, axs = plot_target(target)
print(target.velocity)
print("Import successful, you may proceed.")

domain = "wavelength"       # "frequency" or "wavelength", "wavelength" is recommended
res_type = "log"            # "log" or 'linear', "log" is recommended.
pmin = min(target.wavelength)                    # Minimum value after resampling in units of domain
pmax = max(target.wavelength)                  # Maximum value after resampling in units of domain
pn = 20                     # Number of samples, 20-30 points are recommended.


target.easy_resample(pmin=pmin, pmax=pmax, pn=pn, res_type=res_type, domain=domain, inplace=True)

fig, axs = plot_target(target)
print(target.velocity)

print("Resample successful, you may proceed.")

target_name = "tar5"        # Name of target file (no .target suffix)
target_name=name.split('/')[-1][:-4]
version = "3.4.2"           # Version of Geopsy "2.10.1" or "3.4.2"


# Save to Disk
if os.path.isdir("0_targets/")==False:
    os.mkdir("0_targets/")
target.to_target(f"0_targets/{target_name}", version=version)

# Confirm file exists.
if os.path.exists(f"0_targets/{target_name}.target"):
    print(f"{target_name}.target exists, you may proceed.")

# Minimum and maximum for all parameters. Refer to detailed instructions above.
vp_min, vp_max, vp_dec = 150., 1800., True
vs_min, vs_max, vs_dec = 80., 1000., True
pr_min, pr_max = 0.2, 0.5 ## poisson variable
rh_min, rh_max = 2000., 2000. ## densidad fija

# Layering by Number (LN) parameterizations to consider. Add or remove as desired.
# See Vantassel and Cox (2021) for details.
lns = [3, 4, 5]

# Layering Ratios (LRs) parameterizations to consider. Add or remove as desired.
# See Vantassel and Cox (2021) and Cox and Teague (2016) for details.
lrs = [3.0, 2.0, 1.5]

# Depth factor, typically 2 or 3.
depth_factor = 2

# Minimum and maximum wavelength, selected from experimental disperison data by default.
wmin, wmax = min(target.wavelength), max(target.wavelength)

# Mass density.
if (rh_min - rh_max) < 1:
    rh = swprepost.Parameter.from_fx(rh_min)
else:
    rh = swprepost.Parameter.from_ln(wmin=wmin, wmax=wmax, nlayers=1, par_min=rh_min, par_max=rh_max, par_rev=False)

# Poisson's ratio
if (pr_max - pr_min) < 0.05:
    raise ValueError(
        f"Difference between pr_min and pr_max is too small ({pr_max - pr_min:2f}<0.05), use larger range.")
else:
    pr = swprepost.Parameter.from_ln(wmin=wmin, wmax=wmax, nlayers=1, par_min=pr_min, par_max=pr_max, par_rev=False)

# Make 1_parameters directory.
if not os.path.isdir("1_parameters/"):
    os.mkdir("1_parameters/")

# Parameterize Vs using Layering by Number (LN)
for ln in lns:
    #par_rev dice si el perfil permite inversiones o no
    vs = swprepost.Parameter.from_ln(wmin=wmin, wmax=wmax, nlayers=ln, par_min=vs_min, par_max=vs_max, par_rev=vs_dec,
                                     depth_factor=depth_factor)
    vp = swprepost.Parameter.from_parameter_and_link(par_min=vp_min, par_max=vp_max, par_rev=vp_dec,
                                                     existing_parameter=vs, ptype="vs")
    par = swprepost.Parameterization(vp=vp, pr=pr, vs=vs, rh=rh)
    par.to_param(f"1_parameters/ln{ln}", version=version)

#Parameterize Vs using Layering Ratio (LR)
for lr in lrs:
    vs = swprepost.Parameter.from_lr(wmin=wmin, wmax=wmax, lr=lr, par_min=vs_min, par_max=vs_max, par_rev=vs_dec,
                                     depth_factor=depth_factor)
    vp = swprepost.Parameter.from_parameter_and_link(par_min=vp_min, par_max=vp_max, par_rev=vp_dec,
                                                     existing_parameter=vs, ptype="vs")
    par = swprepost.Parameterization(vp=vp, pr=pr, vs=vs, rh=rh)
    par.to_param(f"1_parameters/lr{int(lr * 10)}", version=version)

nparam = len(lns) + len(lrs)
if len(glob.glob("1_parameters/*.param")) == nparam:
    print(f"All {nparam} .param files exist, you may proceed.")

# Inversion Setup
# ---------------

# Analysis name that is brief, memorable, and descriptive.
# Each output file will begin with this string of characters.
# No spaces or special characters are permitted.
analysis_name =target_name

# Number (positive integer) of inversion trials to perform
# per parameterization. (3 is recommended)
number_of_inversion_trials = 3

# Number (positive integer) of Neighborhood-Algorithm iterations
# to perform per inversion. (250 is recommended)
number_of_iterations = 250

# Number (positive integer) of randomly sampled profiles to attempt
# before the first Neighborhood-Algorithm iteration. (10000 is recommended)
number_of_initial_random_samples = 10000

# Number (positive integer) of best profiles to consider when
# resampling. (100 is recommended)
number_of_profiles_to_consider_when_resampling = 100

# Number (positive integer) of new profiles to consider per
# Neighborhood-Algorithm iteration. (200 is recommended)
number_of_profiles_per_iteration = 200

# Results to Export
# -----------------

# Number of ground models/dispersion curves/ellipticity curves to export
number_of_models_to_export = 100

# Number (positive integer) of Rayleigh and Love wave modes to export.
# If no dispersion curves are desired set both the number of Rayleigh and
# Love modes to 0. (1 is recommended)
number_of_rayleigh_modes_to_export = 1
number_of_love_modes_to_export = 0

# Number (positive float) for minimum amd maximum frequency of exported
# dispersion curve(s) in Hz. Selecting a value slightly less than the
# minimum frequency and a value slighlty greater than the maximum frequency
# of your experimental dispersion data is recommended.
minimum_dispersion_frequency = 1.
maximum_dispersion_frequency = 60.

# Number (positive integer) of frequency points in the exported dispersion
# curve(s). (30 is recommended)
number_of_dispersion_frequency_points = 30

# Number (positive integer) of Rayleigh modes to include in exported ellipticity.
# If no ellipticity curves are desired set this value to 0. (1 is recommended)
number_of_rayleigh_ellipticity_modes_to_export = 0

# Number (positive float) for minimum amd maximum frequency of exported
# Rayleigh wave ellipticity curve(s) in Hz. Selecting a value less than and
# greater than the site's resonant frequency is recommended.
minimum_ellipticity_frequency = 0.2
maximum_ellipticity_frequency = 20.

# Number (positive integer) of frequency points in exported Rayleigh wave
# ellipticity curve(s). (64 is recommended)
number_of_ellipticity_frequency_points = 64

# Job Details
# ---------------

# A recognizable name for this job.
# Name is used solely by DesignSafe-CI & AGAVE/TAPIS.
job_name = target_name

# Maximum job runtime in (HH:MM:SS) format.
# If this time is exceeded the job will be canceled by the job scheduler.
# Each queue has its own associated maximum time, typically 48 hours.
# See Stampede2 (if using Geopsy v2.10.1) or Frontera (if using Geopsy 3.4.2) for queue-specfic details.
runtime = "24:00:00"

full_path='/home/doctor/Doctor/Magister/Tesis/databases/process_data/sw_test_seba/Results_00/Inversion00/3_text/'
ndc = 100  # Number of dispersion curves to import, may use "all".
nrayleigh = 1  # Number of Rayleigh modes to import, may use "all".
nlove = 0  # Number of love modes to import, may use "all".
ngm = 100 # Number of ground models to import, may use "all".

#full_path = "./3_text/"
#full_path = "/home/jupyter/MyData/archive/jobs/job-72721be3-5b8c-40bf-8d75-5dbea8952a23-007/Inversion00/3_text/"
# print(os.listdir(full_path))
fnames = glob.glob(full_path + "*_[dD][cC].txt")
# print(fnames)
fnames = [fname[len(full_path):] for fname in fnames]
fnames.sort(key=lambda x: int(re.findall(r".*l[nr](\d+)_tr?\d+_dc.txt", x.lower())[0]))
# # print(fnames)
## dispersion curves and ground models
dcs, gms = {}, {}
for fname in fnames:
    partype, parnumber, seed = re.findall(r".*(l[nr])(\d+)_tr?(\d+)_dc.txt$", fname.lower())[0]
    #partype: ln or lr
    #parnumber: how many layers or LRs
    #seed: number of the trial, starting from 0
    print(partype,parnumber,seed)
    fname = full_path + fname
#
    # Divide LR by 10
    if partype in ['lr']:
        parnumber = str(int(parnumber) / 10)

    # Save by parameterization
    if partype not in dcs.keys():
        dcs.update({partype: {}})
        gms.update({partype: {}})
        firstpass = True

    # Save by parameterization number
    if parnumber not in dcs[partype].keys():
        dcs[partype].update({parnumber: {}})
        gms[partype].update({parnumber: {}})

    # Save by trial
    print('size',os.path.getsize(fname))
    if os.path.getsize(fname) == 0:
        print(f"fname = {fname}, is empty skipping!")
    else:
        dcs[partype][parnumber].update({seed: swprepost.DispersionSuite.from_geopsy(fname=fname, nsets=ndc,
                                                                                    nrayleigh=nrayleigh, nlove=nlove)})
        try:
            gms[partype][parnumber].update(
                {seed: swprepost.GroundModelSuite.from_geopsy(fname=fname[:-6] + "GM.txt", nmodels=ngm)})
        except FileNotFoundError:
            gms[partype][parnumber].update(
                {seed: swprepost.GroundModelSuite.from_geopsy(fname=fname[:-6] + "gm.txt", nmodels=ngm)})

ncols = len(list(dcs.keys()))
fig, axs = plt.subplots(nrows=1, ncols=ncols, sharey=True, figsize=(3 * ncols, 3), dpi=150)
axs = [axs] if type(axs) != np.ndarray else axs
bestseed = {}
blabel = "Each Trial"
fiter = True
for ax, partype in zip(axs, dcs):
    bestseed.update({partype: {}})
    for parnumber in dcs[partype]:
        seeds, misfits = [], []
        for seed in dcs[partype][parnumber].keys():
            seeds.append(seed)
            ## elige el primer misfit ya que es el más pequeño de todos
            misfits.append(dcs[partype][parnumber][seed].misfits[0])
            ## plotea el ultimo por conveniencia (de todos los trials), pero en realidad se puede plotear cualquiera.
            ax.plot(parnumber, misfits[-1], 'bo', label=blabel, alpha=0.2)
            blabel = None
        bestseed[partype].update({parnumber: seeds[misfits.index(min(misfits))]})
    if fiter:
        fiter = False
        ax.legend()
    ax.set_title("Parameterization Type: " + partype)
axs[0].set_ylabel("Dispersion Misfit, " + "$m_{dc}$")
plt.show()
colors = ["tomato", "orange", "gold", "lightgreen", "skyblue", "cyan", "indigo", "violet"]
ndc = 100  # Number of dispersion curves to plot, may use "all".
nray = 1  # Number of Rayleigh-wave modes to plot, may use "all".
nlov = 0  # Number of Love-wave modes to plot, may use "all".

fig, axs = plt.subplots(ncols=2, sharey=True, figsize=(6, 3), dpi=150)

# Plot the Theoretical Modes of Inversion Ground Models.
color_id = 0
lvstack = None
fvstack = None
for partype in dcs:
    for parnumber in dcs[partype]:
        best = bestseed[partype][parnumber]
        suite = dcs[partype][parnumber][best]
        label = f"{partype}={parnumber} {suite.misfit_repr(nmodels=ndc)}"

        color = colors[color_id]
        for dc_count, dcset in enumerate(suite):
            for mode in range(nray):
                try:
                    dc = dcset.rayleigh[mode]
                    # print(label)
                    axs[1].plot(dc.wavelength, dc.velocity, color=color, label=label, linewidth=0.7)
                    label = None
                    axs[0].plot(dc.frequency, dc.velocity, color=color, label=label, linewidth=0.7)
                except KeyError:
                    print(f"Could not find mode {mode}.")
                #curva_lv = np.vstack((dc.wavelength, dc.velocity)).T
                #curva_fv = np.vstack((dc.frequency, dc.velocity)).T
                #print(curva_lv.shape, curva_fv.shape)
            #try:
                #lvstack = np.vstack((lvstack, curva_lv))
                #fvstack = np.vstack((fvstack, curva_fv))

            #except:
                #lvstack = curva_lv
                #fvstack = curva_fv
            if dc_count + 1 == ndc:
                break
        color_id += 1

# Plot the Experimental Dispersion Curve
ax = axs[0]
target.plot(ax=ax)

ax = axs[1]
target.plot(ax=ax, x="wavelength")
#ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax.set_ylabel("")
#print(lvstack.shape, fvstack.shape)
#filesave = name_cdd.split('.')[0]

#target.to_csv(filesave + '_exp_dispersion.txt')
#np.savetxt(filesave + '_fv_dispersion.txt', lvstack)
#np.savetxt(filesave + '_lv_dispersion.txt', fvstack)

print(dir(dc))
plt.show()

gm = 100  # Number of Vs profiles to plot and consider for Vs uncertainty (see next).
plot_depth = 50  # Maximum plot depth in meters.

fig, ax = plt.subplots(figsize=(2, 4), dpi=150)
color_id = 0
all_gm = []
for partype in gms:
    for parnumber in gms[partype]:
        best = bestseed[partype][parnumber]
        suite = gms[partype][parnumber][best]

        label = f"{partype}={parnumber} {suite.misfit_repr(nmodels=ngm)}"
        for gm in suite[:ngm]:
            all_gm.append(gm)
            ax.plot(gm.vs2, gm.depth, color=colors[color_id], label=label, linewidth=0.7)
            label = None
            #depth_vs = np.vstack((gm.depth, gm.vs2)).T

            # try:
            #     stackstack = np.vstack((stackstack, depth_vs))
            # except:
            #     stackstack = depth_vs
        color_id += 1
    ax.set_ylim(plot_depth, 0)
    ax.set_xlabel('Shear Wave Velocity, Vs (m/s)')
    ax.set_ylabel('Depth (m)')
    #ax.legend(bbox_to_anchor=(1, 0.5), loc='center left')
plt.show()

fig, ax = plt.subplots(figsize=(2, 4), dpi=150)
color_id = 0
all_gm_suite = swprepost.GroundModelSuite.from_list(all_gm)
ddepth, dsigmaln = all_gm_suite.sigma_ln()
ax.plot(dsigmaln, ddepth, linewidth=2)
ax.set_ylim(plot_depth, 0)
ax.set_xlabel(r"$\sigma_{ln,Vs}$")
ax.set_ylabel("Depth (m)")

#np.savetxt(filesave+'_GM.txt',stackstack)
#np.savetxt(filesave+'_uncertainty.txt',np.vstack((ddepth,dsigmaln)))

plt.show()
# Changes are not required below this line
# # ----------------------------------------
#
#