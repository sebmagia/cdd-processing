import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as mpl
import os
from shapely.geometry import Point, Polygon, LineString, box
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely.figures import SIZE, BLUE, GRAY, RED, YELLOW, BLACK, set_limits
import swprepost
import matplotlib

matplotlib.use('Agg')
#def weights_lnvs(file):
def non_decreasing_velocity_profile(depth, velocity):
    """
    Generate a non-decreasing velocity profile.

    Parameters:
        depth (array): Array of depth values.
        velocity (array): Array of velocity values.

    Returns:
        A tuple of two arrays: the modified depth and velocity arrays where the velocity profile is non-decreasing.
    """
    # Compute the cumulative maximum of the velocity array
    velocity_cummax = np.maximum.accumulate(velocity)

    # Return the modified depth and velocity arrays
    return depth, velocity_cummax
def plot_means(mwm,mwm2):

    pass
def gen_means(vels, sigmas, depths,  weights, indexes):
    ii = 0
    vels = np.array(vels)
    print(weights,indexes)
    for w, idx in zip(weights, indexes):
        print(len(idx), ii)
        if len(w) > 1:

            idx = np.array(idx)
            models = vels[idx]
            w = np.array(w).reshape((len(w), 1))
            weighted_model = np.sum(w * models, axis=0)
            sigmas_chosen = [sigmas[i] for i in idx]
            #sigma_weighted=((np.squeeze(np.array(sigmas_chosen)).T)@w)
            maxlen = np.max([len(x) for x in sigmas_chosen])
            # create an empty list to store the padded signals
            maxlen = max(len(sig) for sig in sigmas_chosen)  # find the maximum length of all signals

            # np.pad(sigmas_chosen[1],maxlen-len(sigmas_chosen[1]),mode='constant')
            arrs = np.array(
                [np.pad(sig, (0, maxlen - len(sig)), mode='constant', constant_values=1e+16) for sig in
                 sigmas_chosen]).T
            arrs = np.squeeze(arrs)
            weights_sigma = 1 / arrs
            weights_sum = np.sum(weights_sigma, axis=1)
            normalized_weights = weights_sigma / weights_sum[:, np.newaxis]

            auxiliar = [x * normalized_weights * models.T for x in w]
            weighted_model_2 = np.sum(np.sum(auxiliar, axis=0), axis=1)
            weighted_model_2 = non_decreasing_velocity_profile(depths, weighted_model_2)[1]



            try:
                weighted_models = np.vstack((weighted_models, weighted_model))
                mean_weighted_models = np.vstack((mean_weighted_models, np.mean(weighted_model, axis=-1)))

                weighted_models_2 = np.vstack((weighted_models_2, weighted_model_2))
                mean_weighted_models_2 = np.vstack((mean_weighted_models_2, np.mean(weighted_model_2, axis=-1)))

                #sigmas_weighted = np.vstack((sigmas_weighted, sigma_weighted))
                #mean_sigmas_weighted = np.vstack((mean_sigmas_weighted, np.mean(sigmas_weighted, axis=-1)))

            except ValueError:
                weighted_models = weighted_model
                mean_weighted_models = np.mean(weighted_model, axis=-1)

                weighted_models_2 = weighted_model_2
                mean_weighted_models_2 = np.mean(weighted_model_2, axis=-1)

                #sigmas_weighted = sigma_weighted
                #mean_sigmas_weighted =  np.mean(sigmas_weighted, axis=-1)

            # ax.plot(weighted_model,depth_finer,'r')


        if len(w) == 1:
            print(len(w))
            idx = np.array(idx)
            weighted_model = vels[idx].reshape(vels[idx].shape[1], )
            try:
                weighted_models = np.vstack((weighted_models, weighted_model))
                mean_weighted_models = np.vstack((mean_weighted_models, np.mean(weighted_model, axis=-1)))

                weighted_models_2 = np.vstack((weighted_models_2, weighted_model))
                mean_weighted_models_2 = np.vstack((mean_weighted_models_2, np.mean(weighted_model, axis=-1)))

                #sigmas_weighted = np.vstack((sigmas_weighted, sigma_weighted))
                #mean_sigmas_weighted = np.vstack((mean_sigmas_weighted, np.mean(sigmas_weighted, axis=-1)))

            except ValueError:
                print('xd', weighted_model_2.shape)
                print('xb', sigma_weighted.shape)
                weighted_models = weighted_model
                mean_weighted_models = np.mean(weighted_model, axis=-1)

                weighted_models_2 = weighted_model
                mean_weighted_models_2 = np.mean(weighted_model, axis=-1)

                #sigmas_weighted = np.vstack((sigmas_weighted, sigma_weighted))
                #mean_sigmas_weighted =  np.mean(sigmas_weighted, axis=-1)

                print('duras', sigmas_weighted.shape)
                print('drogas',np.mean(sigmas_weighted, axis=-1).shape)

        if len(w) == 0:
            print('duna')
            weighted_model = np.zeros(len(depths))
            print(weighted_model)
            weighted_model_2 = np.zeros(len(depths))
            sigma_weighted = np.zeros(len(depths))
            try:
                weighted_models = np.vstack((weighted_models, weighted_model))
                mean_weighted_models = np.vstack((mean_weighted_models, -999))

                weighted_models_2 = np.vstack((weighted_models_2, weighted_model_2))
                mean_weighted_models_2 = np.vstack((mean_weighted_models_2, -999))


                #sigmas_weighted = np.vstack((sigmas_weighted, sigma_weighted))
                #print(mean_sigmas_weighted)
                #print(mean_weighted_models)
                #mean_sigmas_weighted = np.vstack((mean_sigmas_weighted, 20.0))
            #except ValueError:
            except UnboundLocalError:
                print('sigwey', sigma_weighted)
                weighted_models = weighted_model
                mean_weighted_models = -999

                weighted_models_2 = weighted_model_2
                mean_weighted_models_2 = -999

                #sigmas_weighted = sigma_weighted
                #mean_sigmas_weighted =  20.0
        ii += 1
    return weighted_models, weighted_models_2, mean_weighted_models, mean_weighted_models_2

def plot_models(vels, sigmas, depths, weights, indexes):
    for w,poly,idx in zip(weights,polygons,indexes):
        print(ii,len(idx))
        #print(gm)
        if len(idx)>1:
            idx=np.array(idx)
            models=vels[idx]
            w = np.array(w).reshape((len(w), 1))
            weighted_model = np.sum(w * models, axis=0)
            sigmas_chosen=[sigmas[i] for i in idx]
            maxlen=np.max([len(x) for x in sigmas_chosen])
             # create an empty list to store the padded signals
            maxlen = max(len(sig) for sig in sigmas_chosen)  # find the maximum length of all signals

            #np.pad(sigmas_chosen[1],maxlen-len(sigmas_chosen[1]),mode='constant')
            arrs=np.array([ np.pad(sig,(0,maxlen-len(sig)),mode='constant',constant_values=1e+16) for sig in sigmas_chosen]).T
            arrs=np.squeeze(arrs)
            #idx=indexes[12]
            weights_sigma=1/arrs
            weights_sum=np.sum(weights_sigma,axis=-1)
            #normalized_weights = weights_sigma / weights_sum
            normalized_weights = weights_sigma / weights_sum[:, np.newaxis]

            auxiliar=[x * normalized_weights * models.T for x in w]
            weighted_model_2=np.sum(np.sum(auxiliar, axis=0), axis=1)
            nondv=non_decreasing_velocity_profile(depths, weighted_model_2)
            #weighted_model_2=np.sum(w[0] * w2 * models.T + w[1] * w2 * models.T,axis=1)

            #weighted_model=
            fig, ax = plt.subplots()
            for mod in models:
                ax.plot(mod, depth_finer, 'r')
                ax.plot(weighted_model, depth_finer, 'k')
                ax.plot(weighted_model_2, depth_finer, 'b')
                ax.plot(nondv[1],depth_finer, 'g')
                ax.set_ylim(50, 0)
                ax.set_xlabel('vs')

            plt.savefig(path_wm + 'rect_'+str(ii).zfill(2)+'.png', format='png', dpi=300, bbox_inches='tight')
            plt.close()

        if len(idx)==1:
            idx = np.array(idx)
            mod = vels[idx].reshape(vels[idx].shape[1],)
            fig, ax = plt.subplots()
            ax.plot(mod, depth_finer, 'r')
            ax.set_ylim(50, 0)
            ax.set_xlabel('vs')
            plt.savefig(path_wm + 'rect_' + str(ii).zfill(2) + '.png', format='png', dpi=300, bbox_inches='tight')
            plt.close()
        ii+=1
def calculate_means(gms, gms_sigma, h1, h2, coords, xmin, xmax, ymin, ymax, wide, length,plotear_models=False,plotear_means=True):
   
    cols, rows, rect_isec, weights, indexes, polygons = create_grid(coords, xmin, xmax, ymin, ymax, wide, length)

    depths = []
    vels = []
    sigmas = []
    for gm, gm_sig in zip(gms, gms_sigma):
        print(gm)
        depth_finer, vel_finer = finer_gm(gm)
        # lim=np.where( (depth_finer<=22) & (depth_finer>=14) )
        lim = np.where((depth_finer <= h2) & (depth_finer >= h1))
        depth_finer = depth_finer[lim]
        vel_finer = vel_finer[lim]
        depths.append(depth_finer)
        vels.append(vel_finer)
        sigma = np.loadtxt(path_gm + gm_sig)
        limsig = np.where((sigma[:, 0] <= h2) & (sigma[:, 0] >= h1))
        sigmas.append(sigma[limsig, 1])
    ## get just one value of depths since all are the same
    depths=depths[0]
    if plotear_models:
        plot_models(vels, sigmas, depths, weights, indexes)

    if plotear_means:
        wm, wm2, mwm, mwm2 =gen_means(vels, sigmas, depths, weights, indexes)
        mwm_rs = mwm.reshape((len(cols) - 1, len(rows) - 1))
        mwm2_rs = mwm2.reshape((len(cols) - 1, len(rows) - 1))
        return mwm2_rs, cols, rows
        # fig, ax = plt.subplots(1,2,figsize=(10,8))
        # im = ax[0].imshow(mwm_rs.T, origin='lower',
        #                 extent=[min(cols), max(cols), min(rows), max(rows)])
        # fig.colorbar(im,ax=ax[0],orientation = 'horizontal', shrink = 0.5, label = '$v_s$ [m/s]')
        # fig.suptitle('pata')
        # plt.show(block=False)
        # for coord in coords:
        #     ax[0].plot([coord[0], coord[2]], [coord[1], coord[3]], 'ro-', linewidth=2.5)
        # #ax[0].colorbar(label='Mean Shear Wave Velocity')
        # ax[0].set_xlabel('X distance [m] ')
        # ax[0].set_ylabel('Y distance [m]')
        # #plt.title('')
        # plt.savefig('Ejempli2.png', format='png', dpi=300, bbox_inches='tight')
        #
        # return wm, wm2, mwm, mwm2
        #plot_means(vels,sigmas,weights,indexes)


def create_grid(coords,xmin,xmax,ymin,ymax,wide,length):
    xx = np.hstack((coords[:, 0], coords[:, 2]))
    #xmin = np.min(xx)
    #xmax = np.max(xx)
    yy = np.hstack((coords[:, 1], coords[:, 3]))
    #ymin = np.min(yy)
    #ymax = np.max(yy)
    print(xmin,xmax,ymin,ymax)
    #wide=5
    #length=5
    cols = list(np.arange(xmin-wide, xmax + wide,wide))
    #cols= np.arange(xmin-wide,xmax+wide,20)
    #rows= np.arange(ymin-length, ymax + length, 10)
    rows = list(np.arange(ymin-length, ymax + length, length))
    fig, ax = plt.subplots()

    polygons=[]
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(Polygon([(x, y), (x + wide, y), (x + wide, y + length), (x, y + length)]))
            plot_polygon(polygons[-1], ax=ax)
    plt.close()

    fig, ax = plt.subplots()
    #for rectangle in polygons:
        #plot_polygon(rectangle,ax=ax,color='y')
    k=0
    #rect_isec = len(polygons)*[[]]
    rect_isec = np.empty((len(polygons), 0)).tolist()
    weights = np.empty((len(polygons), 0)).tolist()
    indexes = np.empty((len(polygons), 0)).tolist()
    lines=[]
    for pair in coords:
        line=LineString([pair[0:2],pair[2:4]])
        lines.append(line)
    for n, line in enumerate(lines):
        print(n)
        plot_line(line, ax=ax, color=BLACK,zorder=10)
        for k,rectangle in enumerate(polygons):
            bool_intersection = line.intersects(rectangle)
            intersection = line.intersection(rectangle)
            if bool_intersection and intersection.geom_type == 'LineString':
                #plot_polygon(rectangle, ax=ax, color='r')

                # inter=line.intersection(rectangle).coords.xy
                # isec_r=np.array([inter[0].tolist(),inter[1].tolist()])
                rect_isec[k].append(intersection.length)
                indexes[k].append(n)
                # p1=isec[0].tolist()
                # p2=isec[1].tolist()

                # rect_isec[k].append()
                print('rectanglex ' + str(k).zfill(2) + ' is intersected in ', intersection, ' by line '+ str(n)  )

            if bool_intersection and intersection.geom_type == 'Point':
                print('rectangleb ' + str(k).zfill(2) + ' is intersected in ', intersection , ' by line '+ str(n) )


                #print(intersection)
    cmap = plt.get_cmap('spring')
    lens=[len(x) for x in rect_isec]
    maxl=max(lens)
    print(lens)
    colors = cmap(np.linspace(0, 1.0, maxl+1))
    print(len(colors))
    #c=0
    print(len(rect_isec),len(lens),len(polygons))
    c=0
    for isec, length, rectangle in zip(rect_isec,lens, polygons):
        print(length)
        plot_polygon(rectangle, ax=ax, color=colors[length] )
        weights[c]=[x/sum(isec) for x in isec]
        c+=1
    plt.close()
    return cols, rows, rect_isec,weights, indexes, polygons

def finer_gm(gm):
        model = swprepost.GroundModel.from_geopsy(fname=path_gm + gm)

        ## round depths to the first decimal
        depth = np.round(model.depth, 1)
        velocity = model.vs2

        ## max depth is 9999 so redefine it to 200
        depth[-1] = 200
        ##obtain unique values of depth, depth increases always so no problem here
        depth_uni = np.unique(depth)

        ## obtain unique values of velocity, velocity may decrease with depth so we make a little fix
        indexes = np.unique(velocity, return_index=True)[1]
        vel_uni = np.array([velocity[index] for index in sorted(indexes)])
        ## define a smaller enough dz for merging vel models
        dz = 0.1
        ##
        depth_finer = np.arange(depth[0], depth[-1] + dz, dz)
        vel_finer = np.zeros(len(depth_finer))
        print(vel_finer)
        print(vel_uni)

        for j in range(len(vel_uni)):
            for i in range(len(depth_finer)):
                if depth_finer[i] < depth_uni[j + 1] and depth_finer[i] >= depth_uni[j]:
                    vel_finer[i] = vel_uni[j]

        return depth_finer, vel_finer

    #pass
path_cdd='/home/doctor/Doctor/Magister/Tesis/databases/process_data/CDD_TOMO_FEB_GS_NEW/INV/'
#coords = np.loadtxt(path_cdd + 'cdd_coordinates.dat')
coords=np.loadtxt(path_cdd+'cdd_coordinates.dat')
coords=np.delete(coords,[8,55],axis=0)
xmin=-5
xmax=30
ymin=-20.5
ymax=25
wide=12 ## largo en x
length=12 ## largo en y
h1 = 0
h2 = 50

path_gm='/home/doctor/Doctor/Magister/Tesis/databases/process_data/GROUND_MODELS_FEB_FIXLAY/'
path_wm='/home/doctor/Doctor/Magister/Tesis/databases/process_data/WEIGHTED_MODELS_FEB_FIXLAY/'
gms=sorted(os.listdir(path_gm))
gms=[x for x in  gms if 'SIGMA' not in x]
gms_sigma=sorted(os.listdir(path_gm))
gms_sigma=[x for x in gms_sigma if 'SIGMA' in x]

h1s=np.arange(0,30,5)
h2s=np.arange(5,35,5)
hkeys=[str(x)+ '-' + str(y) + ' [m]' for (x,y) in zip(h1s,h2s)]
n=0
all_mwm2_rss=[]
szs=[3,4,5,6,7,8,9,10,11,12]
all_mwm2_rss={}
n=0
for h1, h2 in zip(h1s, h2s):
    mwm2_rss = []
    for sz in szs:
        mwm2_rs, cols, rows = calculate_means(gms, gms_sigma, h1, h2, coords, xmin, xmax, ymin, ymax, sz, sz)
        mwm2_rss.append(mwm2_rs)
    all_mwm2_rss[hkeys[n]] = mwm2_rss
    n+=1
key=hkeys[0]
for key in hkeys:
    fig, axs = plt.subplots(5, 2, figsize=(15, 15))
    arry=all_mwm2_rss[key]
    arrf=np.hstack([np.ndarray.flatten(x) for x in arry])
    vmin = np.nanmin(np.abs(arrf))
    vmax = np.nanmax(arrf)
    cc=0
    for i in range(5):
        for j in range(2):
            array = arry[cc]
            array[array < - 900] = np.nan
            ax = axs[i, j]
            im = ax.imshow(array.T, origin='lower', extent=[min(cols), max(cols), min(rows), max(rows)],
                           vmin=vmin, vmax=vmax)
            ax.set_aspect('equal')
            ax.text(0.98, 0.98, key +  ' [m]', horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes, color='black')
            # for coord in coords:
            #     ax.plot([coord[0], coord[2]], [coord[1], coord[3]], 'ro-', linewidth=1.5)
            cc += 1
    cax = fig.add_axes([0.3, 0.05, 0.4, 0.02])
    fig.colorbar(im, cax=cax, orientation='horizontal', label='$v_s$ [m/s]')
    fig.suptitle('Sol' + 'X' + str(wide) + 'Y' + str(length) + '[m]')
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.1, hspace=0.1)
    plt.savefig('Solutions/Sol_NEW' + key + '.png', format='png', dpi=300, bbox_inches='tight')


# #for i in range(len())
# for k in range(len(szs)):
#     fig, axs = plt.subplots(5, 2, figsize=(15, 15))
#     arrys=[all_mwm2_rss[k] for x in all_mwm2_rss]
#     vmin = np.nanmin(np.abs(arry))
#     vmax = np.nanmax(arry)
#     cc=0
#     for i in range(5):
#         for j in range(2):
#             arry=arrys[cc]
#             arry[arry < - 900] = np.nan
#             ax = axs[i,j]
#             im = ax.imshow(arry.T, origin='lower', extent=[min(cols), max(cols), min(rows), max(rows)],
#                            vmin=vmin, vmax=vmax)
#             ax.set_aspect('equal')
#             ax.text(0.98, 0.98, str(int(h1s[c])) + '-' + str(int(h2s[c])) + ' [m]', horizontalalignment='right',
#                     verticalalignment='top',
#                     transform=ax.transAxes, color='black')
#             for coord in coords:
#                 ax.plot([coord[0], coord[2]], [coord[1], coord[3]], 'ro-', linewidth=1.5)
#             cc+=1


    # vmin = np.nanmin(np.abs(np.array(mwm2_rss)))
# vmax = np.nanmax(np.array(mwm2_rss))
#
#
# for i in range(3):
#     for j in range(2):
#         mwm2_rss[c][mwm2_rss[c] < -900] = np.nan
#
#         ax = axs[i,j]
#         im = ax.imshow(mwm2_rss[c].T,origin='lower',extent=[min(cols), max(cols), min(rows), max(rows)], vmin = vmin, vmax = vmax)
#         ax.set_aspect('equal')
#         ax.text(0.98, 0.98, str(int(h1s[c])) + '-' + str(int(h2s[c])) + ' [m]', horizontalalignment='right',
#                 verticalalignment='top',
#                 transform=ax.transAxes, color='black')
#         for coord in coords:
#             ax.plot([coord[0], coord[2]], [coord[1], coord[3]], 'ro-', linewidth=1.5)
#         c+=1
#
# cax = fig.add_axes([0.3, 0.05, 0.4, 0.02])
# fig.colorbar(im, cax=cax, orientation='horizontal', label = '$v_s$ [m/s]')
# fig.suptitle('Sol'+'X'+str(wide)+'Y'+str(length)+ '[m]')
# plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.1, hspace=0.1)
# plt.savefig('Solutions/Sol'+'X'+str(wide)+'Y'+str(length)+'.png', format='png', dpi=300, bbox_inches='tight')
fig, axs = plt.subplots(3, 2, figsize=(15,15))
# c=0
# vmin = np.nanmin(np.abs(np.array(mwm2_rss)))
# vmax = np.nanmax(np.array(mwm2_rss))
#
#
# for i in range(3):
#     for j in range(2):
#         mwm2_rss[c][mwm2_rss[c] < -900] = np.nan
#
#         ax = axs[i,j]
#         im = ax.imshow(mwm2_rss[c].T,origin='lower',extent=[min(cols), max(cols), min(rows), max(rows)], vmin = vmin, vmax = vmax)
#         ax.set_aspect('equal')
#         ax.text(0.98, 0.98, str(int(h1s[c])) + '-' + str(int(h2s[c])) + ' [m]', horizontalalignment='right',
#                 verticalalignment='top',
#                 transform=ax.transAxes, color='black')
#         for coord in coords:
#             ax.plot([coord[0], coord[2]], [coord[1], coord[3]], 'ro-', linewidth=1.5)
#         c+=1
#
# cax = fig.add_axes([0.3, 0.05, 0.4, 0.02])
# fig.colorbar(im, cax=cax, orientation='horizontal', label = '$v_s$ [m/s]')
# fig.suptitle('Sol'+'X'+str(wide)+'Y'+str(length)+ '[m]')
# plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.1, hspace=0.1)
# plt.savefig('Solutions/Sol'+'X'+str(wide)+'Y'+str(length)+'.png', format='png', dpi=300, bbox_inches='tight')


# im = ax[0].imshow(mwm_rs.T, origin='lower',
#                 extent=[min(cols), max(cols), min(rows), max(rows)])
# fig.colorbar(im,ax=ax[0],orientation = 'horizontal', shrink = 0.5, label = '$v_s$ [m/s]')
# fig.suptitle('pata')
# plt.show(block=False)
# for coord in coords:
#     ax[0].plot([coord[0], coord[2]], [coord[1], coord[3]], 'ro-', linewidth=2.5)
# #ax[0].colorbar(label='Mean Shear Wave Velocity')
# ax[0].set_xlabel('X distance [m] ')
# ax[0].set_ylabel('Y distance [m]')
# #plt.title('')
# plt.savefig('Ejempli2.png', format='png', dpi=300, bbox_inches='tight')

#plot_means(vels,sigmas,weights,indexes)# fig, ax = plt.subplots(1,2,figsize=(10,8))
        # im = ax[0].imshow(mwm_rs.T, origin='lower',
        #                 extent=[min(cols), max(cols), min(rows), max(rows)])
        # fig.colorbar(im,ax=ax[0],orientation = 'horizontal', shrink = 0.5, label = '$v_s$ [m/s]')
        # fig.suptitle('pata')
        # plt.show(block=False)
        # for coord in coords:
        #     ax[0].plot([coord[0], coord[2]], [coord[1], coord[3]], 'ro-', linewidth=2.5)
        # #ax[0].colorbar(label='Mean Shear Wave Velocity')
        # ax[0].set_xlabel('X distance [m] ')
        # ax[0].set_ylabel('Y distance [m]')
        # #plt.title('')
        # plt.savefig('Ejempli2.png', format='png', dpi=300, bbox_inches='tight')
        #
        # return wm, wm2, mwm, mwm2
        #plot_means(vels,sigmas,weights,indexes)





# sigma1=np.loadtxt(path_gm+gms_sigma[0])
# sigma2=np.loadtxt(path_gm+gms_sigma[1])
# if len(sigma1) > len(sigma2):
#     sigma2_pad=np.pad(sigma2,((0, len(sigma1)-len(sigma2)), (0, 0)),mode='constant',constant_values=1e16)
#     sigma2_pad[:,0]=sigma1[:,0]
#     sigma12=np.vstack((sigma1[:,1],sigma2_pad[:,1])).T
# if len(sigma1) < len(sigma2):
#     sigma1_pad=np.pad(sigma1,((0, len(sigma2)-len(sigma1)), (0, 0)),mode='constant')
#     sigma1_pad[:,0]=sigma2[:,0]
#     sigma12=np.vstack((sigma1_pad[:,1],sigma2)).T
#
# else:
#     sigma12=np.vstack((sigma1,sigma2)).T
#
# fig,ax=plt.subplots()
#
# models=[]
#
# depths=[]
# vels=[]
# sigmas=[]
# #fig,ax=plt.subplots()
# h1=0
# h2=6
# for gm,gm_sig in zip(gms,gms_sigma):
#     print(gm)
#     depth_finer, vel_finer = finer_gm(gm)
#     #lim=np.where( (depth_finer<=22) & (depth_finer>=14) )
#     lim = np.where((depth_finer<=h2) & (depth_finer>=h1))
#     depth_finer=depth_finer[lim]
#     vel_finer=vel_finer[lim]
#     depths.append(depth_finer)
#     vels.append(vel_finer)
#     sigma=np.loadtxt(path_gm+gm_sig)
#     limsig=np.where((sigma[:,0]<=h2) & (sigma[:,0]>=h1))
#     sigmas.append(sigma[limsig,1])
#     #except Value Error
#     #try:
#     #except:
#      #   sigmas=sigma[:,1]
# #sigmas=sigmas.T
#     #ax.plot(vel_scale,depth_scale,'r')
#     #ax.plot(vel_finer,depth_finer, 'r')
# #ax.set_ylim(25,0)
# #ax.set_xlabel('vs')
# vels=np.array(vels)
# ii=0
# #si=sigmas[:,indexes[12]]
# #weights=1/si
# #weights_sum=np.sum(weights,axis=1)
# #normalized_weights=weights/weights_sum[:,np.newaxis]
# plot_models=False
#
# if plot_models:
#     for w,poly,idx in zip(weights,polygons,indexes):
#         print(ii,len(idx))
#         #print(gm)
#         if len(idx)>1:
#             idx=np.array(idx)
#             models=vels[idx]
#             w = np.array(w).reshape((len(w), 1))
#             weighted_model = np.sum(w * models, axis=0)
#             sigmas_chosen=[sigmas[i] for i in idx]
#             maxlen=np.max([len(x) for x in sigmas_chosen])
#              # create an empty list to store the padded signals
#             maxlen = max(len(sig) for sig in sigmas_chosen)  # find the maximum length of all signals
#
#             #np.pad(sigmas_chosen[1],maxlen-len(sigmas_chosen[1]),mode='constant')
#             arrs=np.array([ np.pad(sig,(0,maxlen-len(sig)),mode='constant',constant_values=1e+16) for sig in sigmas_chosen]).T
#             arrs=np.squeeze(arrs)
#             #idx=indexes[12]
#             weights_sigma=1/arrs
#             weights_sum=np.sum(weights_sigma,axis=-1)
#             #normalized_weights = weights_sigma / weights_sum
#             normalized_weights = weights_sigma / weights_sum[:, np.newaxis]
#
#             auxiliar=[x * normalized_weights * models.T for x in w]
#             weighted_model_2=np.sum(np.sum(auxiliar, axis=0), axis=1)
#             nondv=non_decreasing_velocity_profile(depths, weighted_model_2)
#             #weighted_model_2=np.sum(w[0] * w2 * models.T + w[1] * w2 * models.T,axis=1)
#
#             #weighted_model=
#             fig, ax = plt.subplots()
#             for mod in models:
#                 ax.plot(mod, depth_finer, 'r')
#                 ax.plot(weighted_model, depth_finer, 'k')
#                 ax.plot(weighted_model_2, depth_finer, 'b')
#                 ax.plot(nondv[1],depth_finer, 'g')
#                 ax.set_ylim(50, 0)
#                 ax.set_xlabel('vs')
#
#             plt.savefig(path_wm + 'rect_'+str(ii).zfill(2)+'.png', format='png', dpi=300, bbox_inches='tight')
#             plt.close()
#
#         if len(idx)==1:
#             idx = np.array(idx)
#             mod = vels[idx].reshape(vels[idx].shape[1],)
#             fig, ax = plt.subplots()
#             ax.plot(mod, depth_finer, 'r')
#             ax.set_ylim(50, 0)
#             ax.set_xlabel('vs')
#             plt.savefig(path_wm + 'rect_' + str(ii).zfill(2) + '.png', format='png', dpi=300, bbox_inches='tight')
#             plt.close()
#         ii+=1
#
#
# print('aww')
# #fig,ax=plt.subplots()
# ii=0
# weighted_models=None
# mean_weighted_models=None
# weighted_models_2=None
# mean_weighted_models_2=None
# for w,idx in zip(weights,indexes):
#     print(len(idx),ii)
#
#     if len(w)>1:
#
#         idx = np.array(idx)
#         models = vels[idx]
#         w = np.array(w).reshape((len(w), 1))
#         weighted_model = np.sum(w * models, axis=0)
#         sigmas_chosen = [sigmas[i] for i in idx]
#         maxlen = np.max([len(x) for x in sigmas_chosen])
#         # create an empty list to store the padded signals
#         maxlen = max(len(sig) for sig in sigmas_chosen)  # find the maximum length of all signals
#
#         # np.pad(sigmas_chosen[1],maxlen-len(sigmas_chosen[1]),mode='constant')
#         arrs = np.array(
#             [np.pad(sig, (0, maxlen - len(sig)), mode='constant', constant_values=1e+16) for sig in sigmas_chosen]).T
#         arrs = np.squeeze(arrs)
#         weights_sigma = 1 / arrs
#         weights_sum = np.sum(weights_sigma, axis=1)
#         normalized_weights = weights_sigma / weights_sum[:, np.newaxis]
#
#         auxiliar = [x * normalized_weights * models.T for x in w]
#         weighted_model_2 = np.sum(np.sum(auxiliar, axis=0), axis=1)
#         weighted_model_2=non_decreasing_velocity_profile(depths, weighted_model_2)[1]
#
#         try:
#             weighted_models=np.vstack((weighted_models,weighted_model))
#             mean_weighted_models=np.vstack((mean_weighted_models,np.mean(weighted_model,axis=-1)))
#             weighted_models_2 = np.vstack((weighted_models_2, weighted_model_2))
#             mean_weighted_models_2 = np.vstack((mean_weighted_models_2, np.mean(weighted_model_2, axis=-1)))
#         except ValueError:
#             weighted_models=weighted_model
#             mean_weighted_models=np.mean(weighted_model,axis=-1)
#             weighted_models_2 = weighted_model_2
#             mean_weighted_models_2 = np.mean(weighted_model_2, axis=-1)
#         #ax.plot(weighted_model,depth_finer,'r')
#
#
#     elif len(w)==1:
#         print(len(w))
#         idx = np.array(idx)
#         weighted_model = vels[idx].reshape(vels[idx].shape[1], )
#         try:
#             weighted_models=np.vstack((weighted_models,weighted_model))
#             mean_weighted_models=np.vstack((mean_weighted_models,np.mean(weighted_model,axis=-1)))
#             weighted_models_2 = np.vstack((weighted_models_2, weighted_model))
#             mean_weighted_models_2 = np.vstack((mean_weighted_models_2, np.mean(weighted_model, axis=-1)))
#         except ValueError:
#             weighted_models=weighted_model
#             mean_weighted_models=np.mean(weighted_model,axis=-1)
#             weighted_models_2 = weighted_model
#             mean_weighted_models_2 = np.mean(weighted_model, axis=-1)
#
#     else:
#         weighted_model=np.zeros(len(depth_finer))
#         weighted_model_2=np.zeros(len(depth_finer))
#
#         try:
#             weighted_models = np.vstack((weighted_models, weighted_model))
#             mean_weighted_models = np.vstack((mean_weighted_models, 150.0))
#
#             weighted_models_2 = np.vstack((weighted_models_2, weighted_model_2))
#             mean_weighted_models_2 = np.vstack((mean_weighted_models_2, 150.0))
#         except ValueError:
#             weighted_models = weighted_model
#             mean_weighted_models = 150.0
#
#             weighted_models_2 = weighted_model_2
#             mean_weighted_models_2 = 150.0
#     ii+=1
# #ax.set_ylim(25,0)
# #ax.set_xlabel('vs')
# #plt.show(block=False)
#
# cmap = plt.get_cmap('spring')
# mean_weighted_models_nz= mean_weighted_models[np.nonzero(mean_weighted_models)]
# mean_weighted_models_rs=mean_weighted_models.reshape((len(cols)-1,len(rows)-1))
#
# mean_weighted_models_nz_2= mean_weighted_models_2[np.nonzero(mean_weighted_models_2)]
# mean_weighted_models_rs_2=mean_weighted_models_2.reshape((len(cols)-1,len(rows)-1))
# #norm=(mean_weighted_models-np.min(mean_weighted_models_nz))/(np.max(mean_weighted_models)-np.min(mean_weighted_models_nz))
#
#
#
# fig,ax=plt.subplots()
# im=plt.imshow(mean_weighted_models_rs_2.T,origin='lower',extent=[min(cols),max(cols),min(rows),max(rows)])
# plt.show(block=False)
# for coord in coords:
#     plt.plot([coord[0],coord[2]],[coord[1],coord[3]],'ro-',linewidth=2.5)
# plt.colorbar(label='Mean Shear Wave Velocity')
# plt.xlabel('X distance [m] ')
# plt.ylabel('Y distance [m]')
# plt.savefig('Ejemplito_5.png',format='png', dpi=300, bbox_inches='tight')
#norm = mpl.Normalize(vmin=-1, vmax=1)
#colors = cmap(np.linspace(0, 1.0, maxl+1))
# colors = cmap(norm)
# fig,ax=plt.subplots()
# centroid_x=np.array([x.centroid.x for x in polygons])
# centroid_y=np.array([x.centroid.y for x in polygons])
# center_coords=np.vstack((centroid_x,centroid_y)).T
# for mm,polygon,color in zip(mean_weighted_models,polygons,colors):
#     if mm < 1e-8:
#         plot_polygon(polygon, ax=ax, color='k')
#     else:
#         plot_polygon(polygon,ax=ax,color=color)
# plt.show(block=False)
# # xmin=-5
# # xmax=30
# # ymin=-22
# # ymax=25
# fig,ax=plt.subplots()
# im=plt.imshow(mean_weighted_models_rs.T,origin='lower',extent=[-11,31,-28,26])
# plt.show(block=False)
# for coord in coords:
#     plt.plot([coord[0],coord[2]],[coord[1],coord[3]],'ro-',linewidth=2.5)
# plt.colorbar(label='Mean Shear Wave Velocity')
# plt.xlabel('X distance [m] ')
# plt.ylabel('Y distance [m]')
# plt.savefig('Ejemplito.png',format='png', dpi=300, bbox_inches='tight')

#plt.colorbar(ax=ax,location='left')

#cax=plt.axes([0.85, 0.1, 0.075, 0.8])
#cm=cmap.ScalarMappable(cmap=colors)
#fig.colorbar(cmap.ScalarMappablne(norm=norm, cmap=cmap), ax=ax)    #ax.plot(vel_finer,depth_finer, 'r')
        ## modelos a considerar

    #print(len(w))
#     model=swprepost.GroundModel.from_geopsy(fname=path_gm+gm)
#     depth=np.round(model.depth,2)
#     idx=np.where(depth <= 60.0)
#     depth = depth[idx]
#     models.append(model)
# model=models[0]
# depth=np.round(model.depth,1)
# depth[-1]=200
# depth_uni=np.unique(depth)
# vel_model=model.vs2
#
# depth_scale=np.arange(depth[0],depth[-1]+dz,dz)
# print(depth_scale.shape)
#
# vel_scale=np.zeros(len(depth_scale))


#for i in range(len(depth_scale)):


# fig,ax=plt.subplots()
# ax.plot(vel_scale,depth_scale,'r')
# ax.plot(model.vs2,model.depth,'r')
# ax.set_ylim(60,0)
# ax.set_xlabel('vs')
#
#
# model=models[1]
# depth=np.round(model.depth,1)
# depth[-1]=200
# depth_uni=np.unique(depth)
# vel_model=model.vs2[:-1]
# indexes=np.unique(vel_model,return_index=True)[1]
# vels=[vel_model[index] for index in sorted(indexes)]
# dz=0.1
# depth_scale2=np.arange(depth[0],depth[-1]+dz,dz)
# print(depth_scale.shape)
# vel_scale2=np.zeros(len(depth_scale))
#
#
# #for i in range(len(depth_scale)):
#
# for j in range(len(vels)):
#     for i in range(len(depth_scale2)):
#         if depth_scale2[i] < depth_uni[j+1] and depth_scale[i] >= depth_uni[j]:
#             vel_scale2[i] = vels[j]
#
# #fig,ax=plt.subplots()
# ax.plot(vel_scale2,depth_scale2,'k')
# ax.plot((vel_scale+vel_scale2)/2,depth_scale2,'b')
# ax.plot(model.vs2,model.depth,'k')
# ax.set_ylim(60,0)
# ax.set_xlabel('vs')



# for i in range(len(depth_scale)):
#     for j in range(len(vels)):
#         if depth_scale[i] < depth_uni[j+1]:
#             vel_scale[i]=vels[j]
#idx=np.where(depth <= 60.0)[0]



