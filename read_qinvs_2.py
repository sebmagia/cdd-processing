import numpy as np
import matplotlib.pyplot as plt
import shutil
import swprepost
import os
import matplotlib
matplotlib.use('Agg')



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


    if nlay > 3 and nlay <=6:
        lines=[x.strip() for x in lines if 'generation' not in x]
        lines_1=lines[::2]
        lines_2=lines[1::2]
        lines_11=[" ".join(x.split()) for x in lines_1]
        lines_22=[" ".join(x.split()) for x in lines_2]

        lines_3=[x+" "+y for x,y in zip(lines_11,lines_22)]
        for line in lines_3:
            line = line.split(' ')
            line= np.array([float(x) for x in line])
            #print('line',len(line))
            try:
                models = np.vstack((models, line))
            except NameError:
                models = line

    if nlay ==7:
        lines = [x.strip() for x in lines if 'generation' not in x]
        lines_1 = lines[::3]
        lines_2 = lines[1::3]
        lines_3 = lines[2::3]

        lines_11 = [" ".join(x.split()) for x in lines_1]
        lines_22 = [" ".join(x.split()) for x in lines_2]
        lines_33 = [" ".join(x.split()) for x in lines_3]

        lines_4 = [x + " " + y + " " + z for x, y, z in zip(lines_11, lines_22, lines_33)]
        for line in lines_4:
            line = line.split(' ')
            line = np.array([float(x) for x in line])
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

path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/INV/'
path_GM = '/home/doctor/Doctor/Magister/Tesis/databases/process_data/GROUND_MODELS_FEB_FIXLAY/'

subfolder='qinvs_feb15_fl/'
path_qinv='/home/doctor/Doctor/Magister/Tesis/databases/qinvs/'+subfolder
qinv_folders=sorted(os.listdir(path_qinv))
qinv_folders=[x+'/' for x in qinv_folders if os.path.isdir(path_qinv+x) and x != 'Figures']
cdd_names=os.listdir(path_cdd)
cdd_names=sorted([x for x in cdd_names if '.dat' in x and 'coordinates' not in x])
cdd_names=[x[:-4] for x in cdd_names]

fig_dir = path_qinv + 'Figures/'
for cdd_name in cdd_names[:1]:
    print(cdd_name)

    num = int(cdd_name.split('_')[1])
    cdd=np.loadtxt(path_cdd+cdd_name+'.dat')
    coords=np.loadtxt(path_cdd+'cdd_coordinates.dat')
    folders=sorted([x for x in qinv_folders if cdd_name in x])


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
    all_models=[]
    for folder in folders:
        suite,nlay,seed,uni_models,freq,c_obs,c_cal=obtain_GM(path_qinv+folder,NM=1)
        for gm in suite:
            all_models.append(gm)
        suites.append(suite)
        c_cals.append(c_cal)
        nlays.append(nlay)
        seeds.append(seed)
        uni_modelss.append(uni_models)
    from itertools import chain
    depths=list(chain(*[x.depth for x in all_models]))
    depths=[x for x in depths if x < 7000]
    factor=1.1
    maxdepth=factor*np.max(depths)
    fig,ax=plt.subplots(3,2,figsize=(20,20))
    all_gm_suite = swprepost.GroundModelSuite.from_list(all_models)
    #ax[0,0].set_ylim(50, 0)
    ax[0,0].set_ylim(maxdepth,0)
    ax[0,0].set_xlabel("Vs (m/s)")
    ax[0,0].set_ylabel('Depth [m]')
    #ax[0,1].set_ylim(50, 0)
    ax[0,0].set_ylim(maxdepth,0)
    ax[0,1].set_xlabel("Vs (m/s)")
    ax[0,1].set_ylabel('Depth [m]')
    ax[0,1].set_ylim(maxdepth,0)

    ax[1,1].set_ylabel('$C_R [m/s]$')
    ax[1,1].set_xlabel('Frequency [Hz]')
    ax[1,1].plot(freq,c_obs,'k',linewidth=2,zorder=3)

    ax[1,0].set_ylabel('Misfit [m/s]')
    ax[1,0].set_xlabel('Frequency [Hz]')

    ax[2,1].set_ylabel('Misfit [m/s]')
    ax[2,1].set_xlabel('Number of Layers')
    ddepth, dsigmaln = all_gm_suite.sigma_ln(dmax=np.round(maxdepth),dy=0.1)
    np.savetxt(path_GM+'GM_SIGMA'+cdd_name+'.txt',np.vstack((ddepth,dsigmaln)).T)
    ddepth=np.round(ddepth,2)
    ax[2,0].plot(dsigmaln, ddepth, linewidth=0.75)
    #ax[2,0].set_ylim(50,0)
    ax[2,0].set_ylim(maxdepth,0)
    ax[2,0].set_xlabel("$\sigma_{ln,Vs}$")
    ax[2,0].set_ylabel("Depth[m]")

    max_misfit=max(all_gm_suite.misfits)
    min_misfit=min(all_gm_suite.misfits)

    #from matplotlib import cm
    for cteo,gm in zip(c_cals,all_gm_suite):
        label = 'N' + str(len(gm.vs))
        if len(gm.vs)==7:
            ax[0,0].plot(gm.vs2, gm.depth, color=plt.cm.jet(gm.misfit/(max_misfit) ))
            ax[0,1].plot(gm.vs2, gm.depth, color='green')
            ax[1,1].plot(freq,cteo,color='green')
            ax[1,0].plot(freq,np.abs(cteo-c_obs),color='green')
            ax[2,1].plot(len(gm.vs),gm.misfit,'o',label=label,color='green')

        if len(gm.vs)==6:
            ax[0,0].plot(gm.vs2, gm.depth, color=plt.cm.jet(gm.misfit/(max_misfit) ))
            ax[0,1].plot(gm.vs2, gm.depth, color='tomato')
            ax[1,1].plot(freq,cteo,color='tomato')
            ax[1,0].plot(freq,np.abs(cteo-c_obs),color='tomato')
            ax[2,1].plot(len(gm.vs),gm.misfit,'o',label=label,color='tomato')

        if len(gm.vs)==5:
            ax[0,0].plot(gm.vs2, gm.depth, color=plt.cm.jet(gm.misfit/(max_misfit) ))
            ax[0,1].plot(gm.vs2, gm.depth, color='gold')
            ax[1,1].plot(freq,cteo,color='gold')
            ax[1,0].plot(freq,np.abs(cteo-c_obs),color='gold')
            ax[2,1].plot(len(gm.vs),gm.misfit,'o',label=label,color='gold')

        if len(gm.vs)==4:
            ax[0,0].plot(gm.vs2, gm.depth, color=plt.cm.jet(gm.misfit/(max_misfit) ))
            ax[0,1].plot(gm.vs2, gm.depth, color='cyan')
            ax[1,1].plot(freq,cteo,color='cyan')
            ax[1,0].plot(freq,np.abs(cteo-c_obs),color='cyan')
            ax[2,1].plot(len(gm.vs),gm.misfit,'o',label=label,color='cyan')

        if len(gm.vs)==3:
            ax[0,0].plot(gm.vs2, gm.depth, color=plt.cm.jet(gm.misfit/(max_misfit) ))
            ax[0,1].plot(gm.vs2, gm.depth, color='blue')
            ax[1,1].plot(freq,cteo,color='blue')
            ax[1,0].plot(freq,np.abs(cteo-c_obs),color='blue')
            ax[2,1].plot(len(gm.vs),gm.misfit,'o',label=label,color='blue')

        print(plt.cm.jet(gm.misfit))
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=min_misfit, vmax=max_misfit))
    fig.colorbar(sm,ax=ax[0,0])
    plt.savefig(fig_dir+'/'+cdd_name+'_allinone.png',format='png',dpi=300,bbox_inches='tight')
    plt.close()
    #fig.colorbar(ax=ax)
    #ax.legend()

    fig, ax = plt.subplots(figsize=(3,6), dpi=150)
    for n,suit in enumerate(suites):
        label='N'+str(nlays[n])+'_'+'S'+str(seeds[n])
        gm=suit[0]
        print(nlays[n])
        if nlays[n] == 7:
                ax.plot(gm.vs2, gm.depth,label=label, color = 'green')
        if nlays[n] == 6:
                ax.plot(gm.vs2, gm.depth,label=label, color = 'tomato')
        if nlays[n] == 5:
                ax.plot(gm.vs2, gm.depth, label=label, color = 'gold')
        if nlays[n] == 4:
                ax.plot(gm.vs2, gm.depth, label=label, color = 'cyan')
        if nlays[n] == 3:
                ax.plot(gm.vs2, gm.depth, label=label, color= 'blue')

    mf=[x.misfits for x in suites]
    best_model=suites[np.argmin(mf)]
    best_model.write_to_txt(path_GM+'GM_'+cdd_name+'.txt')
    #ax.legend()
    #ax.set_ylim(50,0)
    ax.set_ylim(maxdepth, 0)
    ax.set_xlabel("Vs (m/s)")
    plt.savefig(fig_dir+'/'+cdd_name+'_models.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(6,6), dpi=150)
    for n,co in enumerate(c_cals):
        label = 'N' + str(nlays[n]) + '_' + 'S' + str(seeds[n])
        if nlays[n] == 7:
            ax.plot(freq, co, label=label, color='green')
        if nlays[n] == 6:
            ax.plot(freq, co, label=label, color='tomato')
        if nlays[n] == 5:
            ax.plot(freq, co, label=label, color='gold')
        if nlays[n] == 4:
            ax.plot(freq, co, label=label, color='cyan')
        if nlays[n] == 3:
            ax.plot(freq, co, label=label, color='blue')


    ax.plot(freq, c_obs, label='OBS', color='k',linewidth=2)
    #ax.legend()
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Phase Velocity [m/s]')
    plt.savefig(fig_dir+'/'+cdd_name+'_fit.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()


    fig, ax = plt.subplots(figsize=(6,6), dpi=150)
    for n,co in enumerate(c_cals):
        label = 'N' + str(nlays[n]) + '_' + 'S' + str(seeds[n])
        if nlays[n] == 7:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='green')
        if nlays[n] == 6:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='tomato')
        if nlays[n] == 5:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='gold')
        if nlays[n] == 4:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='cyan')
        if nlays[n] == 3:
            ax.plot(freq, np.abs(co - c_obs), label=label, color='blue')


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
        if nlays[n] == 7:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='green')
        if nlays[n] == 6:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='tomato')
        if nlays[n] == 5:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='gold')
        if nlays[n] == 4:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='cyan')
        if nlays[n] == 3:
            ax.plot(len(gm.vp),gm.misfit,'o',label=label,color='blue')

    #ax.legend()
    ax.set_xlabel('N Layers')
    ax.set_ylabel('Misfit')


    plt.savefig(fig_dir+'/'+cdd_name+'_modelcomp.png',format='png', dpi=300, bbox_inches='tight')
    plt.close()



