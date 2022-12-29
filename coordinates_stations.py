## ARR1A, 8 SEPTIEMBRE  2022
import numpy as np
import itertools
import matplotlib.pyplot as plt

sept22=True
dic21=True
mar22=True
alf21=True
dic21b=False
mar22b=False
dic22=True

if mar22b:
    MED = 'M1'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [0.0, 0.0]  ## correcto
    E2 = [5, 0.0]
    E3 = [12, 2]
    E4 = [3, 6]  # correcto
    E5 = [7, 7]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'co-')
        plt.plot(x[:, 0], x[:, 1], 'co')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M2'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [0.0, 0.0]  ## correcto
    E2 = [5.0, 0.0]
    E3 = [12, 2]
    E4 = [14,8]  # correcto
    E5 = [2, 11]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'co-')
        plt.plot(x[:, 0], x[:, 1], 'co')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M3'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [0.0, 0.0]  ## correcto
    E2 = [5.0, 0.0]
    E3 = [12, 2]
    E4 = [9, 12]  # correcto
    E5 = [14, 12]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        plt.plot(x[:, 0], x[:, 1], 'co')
    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M4'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [3, 6]  ## correcto
    E2 = [7, 7]
    E3 = [14, 8]
    E4 = [9, 12]  # correcto
    E5 = [14, 12]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'co-')
        plt.plot(x[:, 0], x[:, 1], 'co')
    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M5'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [3, 16]  ## correcto
    E2 = [7, 17]
    E3 = [16, 17]
    E4 = [9, 12]  # correcto
    E5 = [14, 12]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'co-')
        plt.plot(x[:, 0], x[:, 1], 'co')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)


    MED = 'M6'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [3, 6]  ## correcto
    E2 = [7, 7]
    E3 = [14, 8]
    E4 = [9, 12]  # correcto
    E5 = [14, 12]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'co-')
        plt.plot(x[:, 0], x[:, 1], 'co')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M7'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [3, 16]  ## correcto
    E2 = [7, 17]
    E3 = [16, 17]
    E4 = [9, 12]  # correcto
    E5 = [14, 12]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'co-')
        plt.plot(x[:, 0], x[:, 1], 'co')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M8'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2B/' + MED + '/'

    E1 = [2, 11]  ## correcto
    E2 = [9, 12]
    E3 = [12, 2]
    E4 = [3, 6]  # correcto
    E5 = [7, 7]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([4, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'co-')
        plt.plot(x[:, 0], x[:, 1], 'co')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

if dic21b:
    MED = 'M1'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1B/' + MED + '/'

    E1 = [0.0, 0.0]  ## correcto
    E2 = [11.9, 0.0]
    E3 = [23.8, 0.0]
    E4 = [0.0, 11.34]  # correcto
    E5 = [11.9, 11.34]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([0, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'yo-')
        plt.plot(x[:, 0], x[:, 1], 'yo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M2'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1B/' + MED + '/'

    E1 = [0.0, 0.0]  ## correcto
    E2 = [11.9, 0.0]
    E3 = [0.0, 21.34]
    E4 = [11.9, 21.34]  # correcto
    E5 = [23.8, 21.34]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([0, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'yo-')
        plt.plot(x[:, 0], x[:, 1], 'yo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M3'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1B/' + MED + '/'

    E1 = [0.0, 11.34]  ## correcto
    E2 = [11.9, 11.34]
    E3 = [0.0, 21.34]
    E4 = [11.9, 21.34]  # correcto
    E5 = [23.8, 21.34]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([0, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'yo-')
        plt.plot(x[:, 0], x[:, 1], 'yo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M4'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1B/' + MED + '/'

    E1 = [0.0, 11.34]  ## correcto
    E2 = [11.9, 11.34]
    E3 = [23.8, 11.34]
    E4 = [0.0, 0.0]  # correcto
    E5 = [11.9, 0.0]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([0, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'yo-')
        plt.plot(x[:, 0], x[:, 1], 'yo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M5'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1B/' + MED + '/'

    E1 = [23.8, 21.34]  ## correcto
    E2 = [23.8, 0]
    E3 = [23.8, 11.34]
    E4 = [11.9, 21.34]  # correcto
    E5 = [0.0, 21.34]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords + np.array([0, 30])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'yo-')
        plt.plot(x[:, 0], x[:, 1], 'yo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)
if dic21:
    MED='M1B'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/'+MED+'/'

    E1=[0.0,0.0] ## correcto
    E2=[12.8,0.0]
    E3=[25.6,0.0]
    E4=[0.45,11.92] # correcto
    E5=[12.8,11.92]
    coords=np.vstack((E1,E2,E3,E4,E5))
    coords=coords-np.array([0,12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        #plt.plot(x[:,0],x[:,1],'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED='M2A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/'+MED+'/'

    E1=[0.0,0.0] ## correcto
    E2=[12.8,0.0]
    E3=[25.6,0.0]
    E4=[12.8,33.22] # correcto
    E5=[23.25,35.76]
    coords=np.vstack((E1,E2,E3,E4,E5))
    coords=coords-np.array([0,12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        #plt.plot(x[:,0],x[:,1],'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path+'coords_'+MED+'.txt',coords)


    MED = 'M3B'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/' + MED + '/'

    E1 = [24.15, 11.92]  ## correcto
    E2 = [1.45, 23.84]
    E3 = [12.8, 24.72]
    E4 = [24.15, 23.84]  # correcto
    E5 = [0.45, 35.76]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords - np.array([0, 12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED='M4'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/'+MED+'/'

    E1=[24.15,11.92] ## correcto
    E2=[12.8,0.0]
    E3=[0.0,0.0]
    E4=[23.25,35.76] # correcto
    E5=[0.45,35.76]
    coords=np.vstack((E1,E2,E3,E4,E5))
    coords=coords-np.array([0,12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        #plt.plot(x[:,0],x[:,1],'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED='M6B'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/'+MED+'/'

    E1=[1.45,23.84] ## correcto
    E2=[25.6,0.0]
    E3=[12.8,24.72]
    #E4=[23.25,35.76] # correcto
    E5=[24.15,11.92]
    coords=np.vstack((E1,E2,E3,E5))
    coords=coords-np.array([0,12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        #plt.plot(x[:,0],x[:,1],'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED = 'M6C'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/' + MED + '/'

    E1 = [1.45, 23.84]  ## correcto
    E2 = [25.6, 0.0]
    E3 = [12.8, 24.72]
    E4 = [0.45,11.92] # correcto
    E5 = [24.15, 11.92]
    coords = np.vstack((E1, E2, E3,E4, E5))
    coords = coords - np.array([0, 12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M7A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/' + MED + '/'

    E1 = [12.8, 11.92]  ## correcto
    E2 = [12.8, 33.22]
    E3 = [23.25, 35.76]
    E4 = [0.45,11.92] # correcto
    E5 = [24.15, 11.92]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords - np.array([0, 12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M8A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/' + MED + '/'

    E1 = [12.8, 0]  ## correcto
    E2 = [12.8, 33.22]
    E3 = [25.6, 0]
    E4 = [0.45, 35.76]  # correcto
    E5 = [24.15, 23.84]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords - np.array([0, 12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M9A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/' + MED + '/'

    E1 = [0, 0]  ## correcto
    E2 = [1.45, 23.84]
    E3 = [25.6, 0]
    E4 = [12.8, 24.72]  # correcto
    E5 = [24.15, 23.84]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords - np.array([0, 12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M10A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/' + MED + '/'

    E1 = [0.45, 11.92]  ## correcto
    E2 = [1.45, 23.84]
    E3 = [12.8, 11.92]
    E4 = [12.8, 24.72]  # correcto
    E5 = [24.15, 23.84]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords - np.array([0, 12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M11A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR1A/' + MED + '/'

    E1 = [0.45, 11.92]  ## correcto
    E2 = [0.45, 35.76]
    E3 = [12.8, 11.92]
    E4 = [23.25, 35.76]  # correcto
    E5 = [24.15, 23.84]
    coords = np.vstack((E1, E2, E3, E4, E5))
    coords = coords - np.array([0, 12])
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'ro-')
        plt.plot(x[:, 0], x[:, 1], 'ro')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)


if mar22:
    MED='M1P'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/'+MED+'/'

    coordtxt=np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/coords_SUR2A.txt')

    E1 = coordtxt[0]  ## correcto
    E2 = coordtxt[1]
    E3 = coordtxt[2]
    E4 = coordtxt[3]  # correcto
    E5 = coordtxt[4]
    E6 = coordtxt[5]

    #coords = np.vstack((E1, E2, E3, E4, E5))
    coords = np.vstack((E1, E2, E3, E4, E5, E6))

    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'bo-')
        plt.plot(x[:, 0], x[:, 1], 'bo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)


    MED='M2P'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/'+MED+'/'

    coordtxt=np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/coords_SUR2A.txt')

    E1 = coordtxt[6]  ## correcto
    E2 = coordtxt[7]
    E3 = coordtxt[8]
    E4 = coordtxt[9]  # correcto
    E5 = coordtxt[10]
    E6 = coordtxt[11]

    #coords = np.vstack((E1, E2, E3, E4, E5))
    coords = np.vstack((E1, E2, E3, E4, E5, E6))

    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'bo-')
        plt.plot(x[:, 0], x[:, 1], 'bo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)


    MED='M3P'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/'+MED+'/'

    coordtxt=np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/coords_SUR2A.txt')

    E1 = coordtxt[6]  ## correcto
    E2 = coordtxt[7]
    E3 = coordtxt[8]
    E4 = coordtxt[12]  # correcto
    E5 = coordtxt[13]
    E6 = coordtxt[14]

    #coords = np.vstack((E1, E2, E3, E4, E5))
    coords = np.vstack((E1, E2, E3, E4, E5, E6))

    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'bo-')
        plt.plot(x[:, 0], x[:, 1], 'bo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M4P'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/' + MED + '/'

    coordtxt = np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/coords_SUR2A.txt')

    E1 = coordtxt[6]  ## correcto
    E2 = coordtxt[7]
    E3 = coordtxt[8]
    E4 = coordtxt[0]  # correcto
    E5 = coordtxt[1]
    E6 = coordtxt[2]

    #coords = np.vstack((E1, E2, E3, E4, E5))
    coords = np.vstack((E1, E2, E3, E4, E5, E6))

    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'bo-')
        plt.plot(x[:, 0], x[:, 1], 'bo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M5P'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/' + MED + '/'

    coordtxt = np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/coords_SUR2A.txt')

    E1 = coordtxt[6]  ## correcto
    E2 = coordtxt[7]
    E3 = coordtxt[8]
    E4 = coordtxt[3]  # correcto
    E5 = coordtxt[4]
    E6 = coordtxt[5]

    #coords = np.vstack((E1, E2, E3, E4, E5))
    coords = np.vstack((E1, E2, E3, E4, E5, E6))

    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'bo-')
        plt.plot(x[:, 0], x[:, 1], 'bo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M6P'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/' + MED + '/'

    coordtxt = np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/coords_SUR2A.txt')

    E1 = coordtxt[9]  ## correcto
    E2 = coordtxt[10]
    E3 = coordtxt[11]
    E4 = coordtxt[3]  # correcto
    E5 = coordtxt[4]
    E6 = coordtxt[5]

    #coords = np.vstack((E1, E2, E3, E4, E5))
    coords = np.vstack((E1, E2, E3, E4, E5, E6))

    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'bo-')
        plt.plot(x[:, 0], x[:, 1], 'bo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'M7P'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/' + MED + '/'

    coordtxt = np.loadtxt('/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/SUR2A/coords_SUR2A.txt')

    E1 = coordtxt[9]  ## correcto
    E2 = coordtxt[10]
    E3 = coordtxt[11]
    E4 = coordtxt[12]  # correcto
    E5 = coordtxt[13]
    E6 = coordtxt[14]

    #coords = np.vstack((E1, E2, E3, E4, E5))
    coords = np.vstack((E1, E2, E3, E4, E5, E6))

    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'bo-')
        plt.plot(x[:, 0], x[:, 1], 'bo')

    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    np.savetxt(path + 'coords_' + MED + '.txt', coords)

if sept22:
    MED='ARR1A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[10.48,2.15]
    E2=[0.0,4.30]
    E3=[20.95,4.30]
    E4=[20.95,0.0]
    E5=[0.0,0.0]

    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')


    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED='ARR1B'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[10.48,0.0]
    E2=[5.24,1.08]
    E3=[15.72,1.08]
    E4=[20.95,0.0]
    E5=[0.0,0.0]
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')

    np.savetxt(path+'coords_'+MED+'.txt',coords)


    MED='ARR1C'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[10.48,4.30]
    E2=[5.24,3.22]
    E3=[15.72,3.22]
    E4=[20.95,0.0]
    E5=[0.0,0.0]
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED='ARR2A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[10.48,15.0] # good
    E2=[2.62,3.76]
    E3=[18.33,3.76]
    E4=[20.95,0.0] # good
    E5=[0.0,0.0] # good
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED='ARR2C'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[10.48,15.0]
    E2=[5.24,7.52]
    E3=[15.71,7.52]
    E4=[20.95,0.0]
    E5=[0.0,0.0]
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')

    np.savetxt(path+'coords_'+MED+'.txt',coords)


    MED='ARR2C'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[10.48,15.0]
    E2=[7.86,11.25]
    E3=[13.09,11.25]
    E4=[20.95,0.0]
    E5=[0.0,0.0]
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED='ARR3A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[23.95,-4.0]
    E2=[23.95,0.0]
    E3=[20.95,3.4]
    E4=[20.95,0.0] # correcto
    E5=[23.95,3.4]
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')
    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED='ARR3B'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[23.95,-4.0] ## correcto
    E2=[23.95,10.2]
    E3=[20.95,6.8]
    E4=[20.95,0.0] # correcto
    E5=[23.95,6.8]
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')
    np.savetxt(path+'coords_'+MED+'.txt',coords)


    MED='ARR3C'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[23.95,-4.0] ## correcto
    E2=[20.95,10.2]
    #E3=[20.95,6.8]
    E4=[20.95,0.0] # correcto
    E5=[23.95,13.6]
    coords=np.vstack((E1,E2,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')
    np.savetxt(path+'coords_'+MED+'.txt',coords)




    MED='ARR3D'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/'+MED+'/'

    E1=[23.95,-4.0] ## correcto
    E2=[20.95,17.0]
    E3=[20.95,13.6]
    E4=[20.95,0.0] # correcto
    E5=[23.95,17.0]
    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        plt.plot(x[:,0],x[:,1],'ko')
    np.savetxt(path+'coords_'+MED+'.txt',coords)
    #plt.savefig('Available_correlations.png',format='png', dpi=300)

if alf21:
    MED = 'Lin_1_20'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/Los_Presidentes_HS/' + MED + '/'
    E1=[0.0,0.0] ## correcto
    E2=[0.0,5.0]
    E3=[0.0,10.0]
    E4=[0.0,15.0] # correcto
    E5=[0.0,20.0]

    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        #plt.plot(x[:,0],x[:,1],'mo-')
        plt.plot(x[:,0],x[:,1],'mo')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED = 'Rect_1_27'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/Los_Presidentes_HS/' + MED + '/'
    E1 = [-3.0, 20.0]  ## correcto
    E2 = [-3.0, -5.0]
    E3 = [24.0, -5.0]
    E4 = [24.0, 20.0]  # correcto
    E5 = [10.5, 20.0]

    coords = np.vstack((E1, E2, E3, E4, E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'mo-')
        plt.plot(x[:,0],x[:,1],'mo')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    MED = 'Lin_2_20'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/Los_Presidentes_HS/' + MED + '/'
    E1 = [0.0, -12.0]  ## correcto
    E2 = [5.0, -12.0]
    E3 = [10.0, -12.0]
    E4 = [15.0, -12.0]  # correcto
    E5 = [20.0, -12.0]

    coords = np.vstack((E1, E2, E3, E4, E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'mo-')
        plt.plot(x[:,0],x[:,1],'mo')

    plt.xlim([-20,40])
    #plt.savefig('ALL_COMBINATIONS.png',dpi=300,bbox_inches='tight')
    np.savetxt(path + 'coords_' + MED + '.txt', coords)

if dic22:
    MED = 'MED1A'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/LP_DIC22/' + MED + '/'
    E1=[0.0,-20.2] ## correcto
    E2=[13.7,-20.2]
    E3=[0.0,-12.0]
    E4=[13.7,-12.0] # correcto
    E5=[6.85,-16.1]

    coords=np.vstack((E1,E2,E3,E4,E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x=np.vstack(x)
        #plt.plot(x[:,0],x[:,1],'mo-')
        plt.plot(x[:,0],x[:,1],'co')

    np.savetxt(path+'coords_'+MED+'.txt',coords)

    MED = 'MED1B'
    path = '/home/doctor/Doctor/Magister/Tesis/databases/LosPresidentes_0922/LP_DIC22/' + MED + '/'
    E1 = [27.4, -20.2]  ## correcto
    E2 = [13.6, -20.2]
    E3 = [27.4, -12.0]
    E4 = [13.7, -12.0]  # correcto
    E5 = [20.55, -16.1]

    coords = np.vstack((E1, E2, E3, E4, E5))
    combs = list(itertools.combinations(coords, 2))
    for x in combs:
        x = np.vstack(x)
        #plt.plot(x[:, 0], x[:, 1], 'mo-')
        plt.plot(x[:,0],x[:,1],'co')

    np.savetxt(path + 'coords_' + MED + '.txt', coords)

    plt.savefig('ALL_COMBINATIONS.png',dpi=300,bbox_inches='tight')
    np.savetxt(path + 'coords_' + MED + '.txt', coords)