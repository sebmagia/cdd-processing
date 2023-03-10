import numpy as np
import matplotlib.pyplot as plt
import shutil
import swprepost
import os



#files=[x for x in files if '.sh' not in x]
#fold=folder_output+files[0]+'/'

def obtain_GM(fold,NM=50):
    seed=fold.split('seed')[1][:-1]

    ## read the log to retrieve number of layers
    with open(fold+'log.dat','r') as fileobject:
        lines=fileobject.readlines()
    nlay=0
    for line in lines:
        if 'VS' in line:
          nlay+=1
   ## read mod_gen to get access to CDD output
    with open(fold + 'mod_gen.dat', 'r') as fileobject:
        lines = fileobject.readlines()

    for n,line in enumerate(lines):
        if 'F(Hz)' in line:
            break
    best_sol=lines[n+1:]
    best_sol=[" ".join(x.split()) for x in best_sol]
    ## list for calculated dispersion curves
    for line in best_sol:
        line = line.split(' ')[:5]
        #print(line)
        line =np.array([float(x) for x in line])

        try:
            output = np.vstack((output, line))
        except NameError:
            output = line
    freq = output[:,1]
    c_obs = output[:,2]
    c_cal = output[:,3]

    ## read mod_gen-all to get ground models data
    with open(fold + 'mod_gen-all.dat', 'r') as fileobject:
        lines = fileobject.readlines()

    print(nlay)

    if nlay == 3:
        lines=[x.strip() for x in lines if 'generation' not in x]
        lines = [" ".join(x.split()) for x in lines]
        for line in lines:
            line = line.split(' ')
            line = np.array([float(x) for x in line])
            try:
                models = np.vstack((models, line))
            except NameError:
                models = line
            # if 'generation' not in line:
            #     line = line.strip().split('  ')
            #     line= np.array([float(x) for x in line])
            #     try:
            #         models = np.vstack((models, line))
            #     except NameError:
            #         models = line


    if nlay > 3:
        lines=[x.strip() for x in lines if 'generation' not in x]
        lines_1=lines[::2]
        lines_2=lines[1::2]
        lines_11=[" ".join(x.split()) for x in lines_1]
        lines_22=[" ".join(x.split()) for x in lines_2]

        lines_3=[x+" "+y for x,y in zip(lines_11,lines_22)]
        for line in lines_3:
            line = line.split(' ')
            line= np.array([float(x) for x in line])
            try:
                models = np.vstack((models, line))
            except NameError:
                models = line





    unique_models=np.unique(models[:,0],return_index=True)
    uni_models=models[unique_models[1],:]
    ## get the N best models
    uni_models=uni_models[:NM,:]
    hs_depth=80*np.ones((len(uni_models),1))
    uni_models=np.hstack((uni_models,hs_depth))

    vp=1300
    rho=2000

    for n,model in enumerate(uni_models):
        if n==0:
            vels=model[1::2]
            dz=model[0::2][1:]
            gm=swprepost.GroundModel(thickness=dz, vs=vels,vp=vels*1.7,density=rho*np.ones(len(dz)),misfit=model[0])
            suite = swprepost.GroundModelSuite(groundmodel=gm)
        else:
            vels = model[1::2]
            dz = model[0::2][1:]
            gm = swprepost.GroundModel(thickness=dz, vs=vels, vp=vels*1.7, density=rho * np.ones(len(dz)),misfit=model[0])
            suite.append(groundmodel=gm)

    return suite,nlay, seed, uni_models,freq, c_obs, c_cal

path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO/'
cdd_names=[0,1,3,4,5,6,8,9,12,13,14,15,16,17,19,20,21,22,23,24,30,32,32,35,37,38,41,46]
#cdd_names=[0]
           #20,21,22,23,24,30,31,32,35,37,
#           38,40,41,42,43,45,46,47,48,49,50,51,53,54,55,56,57,58,59,60,61,62,63]
cdd_names=['cdd'+str(x).zfill(2) for x in cdd_names]
#cdd_names=[cdd_names[0]]
subfolder='qinvs_thursday'
for cdd_name in cdd_names:
    print(cdd_name)
    #cdd_name='cdd00'
    num=int(cdd_name.split('dd')[1])
    cdd=np.loadtxt(path_cdd+cdd_name+'.dat')
    coords=np.loadtxt(path_cdd+'cdd_coordinates.dat')
    folder_output='/home/doctor/Doctor/Magister/Tesis/databases/qinvs/' + subfolder + '/' + cdd_name +'/'
    files=sorted(os.listdir(folder_output))
    folders=[folder_output + x + '/'  for x  in files if  x!='script.sh']


    fig_dir=folder_output.split(cdd_name)[0]+'Figures/'
    fig,ax=plt.subplots()
    dist=((coords[num,0]-coords[num,2])**2 + (coords[num,1]-coords[num,3])**2)**0.5
    for x in coords:
        ax.plot([x[0], x[2]], [x[1], x[3]],'k^')
    ax.plot([coords[num,0],coords[num,2]],[coords[num,1],coords[num,3]],'r-^')
    plt.savefig(fig_dir+'/'+cdd_name+'_pos.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()

    colors = ["tomato", "gold", "cyan", "violet"]

    suites=[]
    nlays=[]
    seeds=[]
    uni_modelss=[]
    c_cals=[]
    for folder in folders:
        suite,nlay,seed,uni_models,freq,c_obs,c_cal=obtain_GM(folder,NM=1)
        suites.append(suite)
        c_cals.append(c_cal)
        nlays.append(nlay)
        seeds.append(seed)
        uni_modelss.append(uni_models)

    #cmap = plt.get_cmap('jet')
    #colors = cmap(np.linspace(0, 1.0, len(suites)))
    #colors=['red','magenta']

    fig, ax = plt.subplots(figsize=(3,6), dpi=150)
    for n,suit in enumerate(suites):
        label='N'+str(nlays[n])+'_'+'S'+str(seeds[n])
        gm=suit[0]
        print(nlays[n])
        if nlays[n] == 6:
                ax.plot(gm.vs2, gm.depth,label=label, color = 'tomato')
        if nlays[n] == 5:
                ax.plot(gm.vs2, gm.depth, label=label, color = 'gold')
        if nlays[n] == 4:
                ax.plot(gm.vs2, gm.depth, label=label, color = 'cyan')
        if nlays[n] == 3:
                ax.plot(gm.vs2, gm.depth, label=label, color= 'tomato')

    mf=[x.misfits for x in suites]
    best_model=suites[np.argmin(mf)]
    path_GM='/home/doctor/Doctor/Magister/Tesis/databases/process_data/GROUND_MODELS/'
    best_model.write_to_txt(path_GM+'GM_'+cdd_name+'.txt')
    #ax.legend()
    ax.set_ylim(50,0)
    ax.set_xlabel("Vs (m/s)")
    plt.savefig(fig_dir+'/'+cdd_name+'_models.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(6,6), dpi=150)
    for n,co in enumerate(c_cals):
        label = 'N' + str(nlays[n]) + '_' + 'S' + str(seeds[n])
        if nlays[n] == 6:
            ax.plot(freq, co, label=label, color='tomato')
        if nlays[n] == 5:
            ax.plot(freq, co, label=label, color='gold')
        if nlays[n] == 4:
            ax.plot(freq, co, label=label, color='cyan')
        if nlays[n] == 3:
            ax.plot(freq, co, label=label, color='tomato')


    ax.plot(freq, c_obs, label='OBS', color='k',linewidth=2)
    #ax.legend()
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Phase Velocity [m/s]')
    plt.savefig(fig_dir+'/'+cdd_name+'_fit.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()


    fig, ax = plt.subplots(figsize=(6,6), dpi=150)
    for n,co in enumerate(c_cals):
        label = 'N' + str(nlays[n]) + '_' + 'S' + str(seeds[n])
        if nlays[n] == 6:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='tomato')
        if nlays[n] == 5:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='gold')
        if nlays[n] == 4:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='cyan')
        if nlays[n] == 3:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='tomato')


    #ax.plot(freq, c_obs, label='OBS', color='k',linewidth=2)
    #ax.legend()
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Misfit of Phase Vel. [m/s]')
    plt.savefig(fig_dir+'/'+cdd_name+'_residuals.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(6,6), dpi=150)
    for n,suit in enumerate(suites):
        label='N'+str(nlays[n])+'_'+'S'+str(seeds[n])
        gm=suit[0]
        if nlays[n] == 6:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='tomato')
        if nlays[n] == 5:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='gold')
        if nlays[n] == 4:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='cyan')
        if nlays[n] == 3:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='tomato')

    #ax.legend()
    ax.set_xlabel('N Layers')
    ax.set_ylabel('Misfit')


    plt.savefig(fig_dir+'/'+cdd_name+'_modelcomp.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()



#suite, nlay, seed, uni_models = obtain_GM(fold,NM=10)
#fold=folder_output+files[1]+'/'
#suite2, nlay2, seed2, uni_models2 = obtain_GM(fold,NM=10)


# fig, ax = plt.subplots(figsize=(2,4), dpi=150)
# # Plot 100 best
# label = "100 Best"
# for gm in suite[:100]:
#     ax.plot(gm.vs2, gm.depth, color="#ababab", label=label)
#     label=None
# # Plot the single best in different color
# #median_dz, median_vs = suite.median_simple(nbest=10, parameter="vs")
#
# ax.plot(suite[0].vs2, suite[0].depth, color="#00ffff", label="1 Best")
# #ax.plot(median_vs, median_dz, color="r", label="median")
#
# ax.set_ylim(50,0)
# ax.set_xlabel("Vs (m/s)")
#
# ## plot uncertainty
# fig, ax = plt.subplots(figsize=(2,4), dpi=150)
# disc_depth, siglnvs = suite[:100].sigma_ln()
# ax.plot(siglnvs, disc_depth, color="#00ff00")
# ax.set_xlim(0, 0.2)
# ax.set_ylim(50,0)
# ax.set_xlabel("$\sigma_{ln,Vs}$")
# ax.set_ylabel("Depth (m)")
# plt.show()
#
# ax.set_ylabel("Depth (m)")
# ax.legend()
# plt.show()


# get better N models
#N=10
#best_models=uni_models[:N,:]
#gm1 = swprepost.GroundModel(thickness=[1.0,0], vs=[100,250],vp=[1300,1300],density=[2000,2000])
#suite=swprepost.GroundModelSuite(groundmodel=gm1)
#suite.append(gm1)
#idx=np.argsort(models[:,0])
# get better N models

#cdd_name='cdd01_ff'