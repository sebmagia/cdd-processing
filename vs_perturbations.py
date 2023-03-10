# =============================================================================

# 4. Generate a random field based on soil element midpoint

# =============================================================================


# define the correlation range in meters for vertical directions

aX = r

# define variance for covariance (variance of ln of Vs)

var = sd ** 2

# make the variance drop to zero in the bottom element and in the numBoundaryEles elements around the perimeter

varList = np.ones(len(eleMidPointX)) * var

####Define width of horizontal and vertical exponential taper to zero variance

expFuncConstHor = 20.0

expFuncConstVer = 2.0

# define width of zero variance zone

zeroVarZoneHor = 100.0

# apply the explonential tpaer to left side and right side of model

xLowVar = 1 - np.e ** (-(eleMidPointX - zeroVarZoneHor) / expFuncConstHor)

xHighVar = 1 - np.e ** (-(numEleX * eleWidth - eleMidPointX - zeroVarZoneHor) / expFuncConstHor)

# apply the exponential taper to base of model

zLowVar = 1 - np.e ** (-(eleMidPointZ - 1.25) / expFuncConstVer)

# get the enveloped variance that should be applied to all points

varFactors = np.array([xLowVar, xHighVar, zLowVar])

varFactors = np.min(varFactors, axis=0)

varFactors = np.array([np.max([i, 0]) for i in varFactors])

varList = varFactors * var

paramsFile.write('zero variance zone width = %s\n' % zeroVarZoneHor)

paramsFile.write('exponential variance decay function width  = %s\n' % expFuncConstHor)

# use GSTools package to create random field of Vs perturbations (applied to Vs profile next)

'''Outdated code (GStools updated??)

cov_model = {'dim': 2, 'mean': .0, 'var': 1.0, 'len_scale': aX,\

             'model': 'exp', 'anis': [anis], 'mode_no': 2000}

srf = SRF(**cov_model)

randField = srf(eleMidPointX, eleMidPointZ, mesh_type='unstructured')

'''

# model = Gaussian(dim=2, mean=0, var=1, len_scale=[aX, aX / anis])

# model = Matern(dim=2, mean=0, var=1, len_scale=[aX, aX / anis], nu=0.2)

model = Exponential(dim=2, var=1, len_scale=[aX, aX / anis])

srf = SRF(model, mean=0)

# srf = SRF(model, seed=int(It))

# srf = SRF(model, seed=19970221)

randField = srf([eleMidPointX, eleMidPointZ])

# scale the random field to the desired variance

VsPerturbations = np.sqrt(varList) * randField

VsPert4Plots = np.sqrt(var) * randField

if var < 0.005:
    VsPerturbations = VsPerturbations * 0.0

# plot probability density function of the entire Vs perturbations random field

plt.figure()

sns.kdeplot(VsPerturbations, color='k', label='Boundary Modified')

sns.kdeplot(VsPert4Plots, color='r', ls='--', label='Unmodified')

plt.xlabel('Vs Pertubation')

plt.ylabel('Probability Density')

ax = plt.gca()

ax.text(0.05, 0.9, '$\mu = $%.2f\n$\sigma_{ln} = $%.2f' % (np.mean(VsPert4Plots), np.std(VsPert4Plots)), ha='left',
        va='top', transform=ax.transAxes)

plt.title('Full model random field Vs perturbations: range = %s, sd = %s' % (r, sd))

plt.legend()

plt.tight_layout()

pdf.savefig()

plt.close()

H
# =============================================================================

# 5. Define Soil Materials

# =============================================================================

# create one material for each element based on the pressure dependent velocity profile developed

inputFile.write('#Create soil materials\n')

eleVs = np.zeros(numEleTot)

eleG = np.zeros(numEleTot)

eleE = np.zeros(numEleTot)

eleRho = np.zeros(numEleTot)

eleVsPerturbed = np.zeros(numEleTot)

eleVsPerturbed4Plots = np.zeros(numEleTot)

# create lists of perturbed Vs values for each individual layer (list of numLayer lists)

indLayerPerturbedVs = []

indLayerPerturbedVs4Plots = []

indLayerPerturbations = []

for i in np.arange(0, numLayers + 1):
    indLayerPerturbedVs.append([])

    indLayerPerturbedVs4Plots.append([])

    indLayerPerturbations.append([])

counter = 0

for material in np.arange(0, numEleTot):
    # find the Vs value from the pressure dependent Vs profile at the depth of each element midpoint

    eleVs[material] = VsDepend[find_nearest(VsZ, [eleMidPoint_Dic[material + 1][1]])][0]

    # perturb the median Vs at that depth (includes zero-variance zone)

    eleVsPerturbed[material] = eleVs[material] * np.exp(VsPerturbations[material])

    # perturb median Vs at that depth (does not include zero variance zones; to calculate distribution parameters without zero variance)

    eleVsPerturbed4Plots[material] = eleVs[material] * np.exp(VsPert4Plots[material])

    indLayerPerturbedVs[eleSoilLayer_Dic[material + 1]].append(eleVsPerturbed[material])

    indLayerPerturbedVs4Plots[eleSoilLayer_Dic[material + 1]].append(eleVsPerturbed4Plots[material])

    eleRho[material] = rho[eleSoilLayer_Dic[material + 1]]

    eleG[material] = eleRho[material] * eleVsPerturbed[material] ** 2.0

    eleE[material] = 2.0 * eleG[material] * (1.0 + nu)

    # create material

    inputFile.write('nDMaterial ElasticIsotropic %s  %s %s %s\n' % (
    material + 1, eleE[material], nu, rho[eleSoilLayer_Dic[material + 1]]))

    # populate Vs Values for GiD file

    gidVsFile.write('%s %s\n' % (material + 1, eleVsPerturbed4Plots[material]))

    gidVsFile_mod.write('%s %s\n' % (material + 1, eleVsPerturbed[material]))

print('mean Vs value = %s' % np.mean(indLayerPerturbedVs4Plots[1]))

print('Std dev ln(Vs) value = %s' % np.std(np.log(indLayerPerturbedVs4Plots[1])))

print('Modified mean Vs value = %s' % np.mean(indLayerPerturbedVs[1]))

print('Modified Std dev ln(Vs) value = %s' % np.std(np.log(indLayerPerturbedVs[1])))