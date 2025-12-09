import numpy as np
import cv2
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from time import time
from scipy.interpolate import interp1d
from IPython.display import clear_output
from scipy.interpolate import UnivariateSpline
import os

x = "x"
y = "y"
z = "z"
absolut = "Absolut"

def vis_acceleration(filnavn, akse):
    
    data = np.loadtxt(filnavn,skiprows=1,delimiter=",")
    tid = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    absolut = data[:,4]
    
    total = 0
    first_second_index = np.where(tid < 1)[0][-1]
    
    if akse == "x":
        for i in x[:first_second_index]:
            total += i
        mean =total/len(x[:first_second_index])
        acceleration = []
        for i in x:
            acceleration.append(i-mean)
    
    if akse == "y":
        for i in y[:first_second_index]:
            total += i
        mean =total/len(y[:first_second_index])
        acceleration = []
        for i in y:
            acceleration.append(i-mean)
    
    if akse == "z":
        for i in z[:first_second_index]:
            total += i
        mean =total/len(z[:first_second_index])
        acceleration = []
        for i in z:
            acceleration.append(i-mean)
    
    if akse == "Absolut":
        for i in absolut[:first_second_index]:
            total += i
        mean =total/len(absolut[:first_second_index])
        acceleration = []
        for i in absolut:
            acceleration.append(i-mean)
    
        
    
    fig, ax = plt.subplots(1,1,figsize = (12,4))
    
    Nyacceleration = acceleration
    ax.plot(tid,Nyacceleration, linewidth = 1, color = "black")
    ax.set_title(r"Acceleration (m/s$^2$)")
    ax.set_xlabel("Tid (s)")
    ax.set_ylabel(r"Acceleration (m/s$^2$)") 
    
    print(f"Maksimal acceleration: {np.max(Nyacceleration):.3f} m/s\u00B2 \n")
    print(f"Maksimal deacceleration: {np.min(Nyacceleration):.3f} m/s\u00B2 \n ")
    
    return tid,acceleration

def vis_bevægelsesligningerne(acceleration):
    
    hastighedsgraf = []
    stedgraf = []
    hastighed_global = 0
    sted_global = 0
    tid_forrige = acceleration[0][0]

    for i,j in zip(acceleration[1],acceleration[0]):
        dt = j-tid_forrige
        hastighed_global+=i*dt
        hastighedsgraf.append(hastighed_global)
        tid_forrige = j

    for i,j in zip(hastighedsgraf,acceleration[0]):
        dt = j-tid_forrige
        sted_global+=i*dt
        stedgraf.append(sted_global)
        tid_forrige = j
    
    fig, ax = plt.subplots(3,1,figsize = (12,12))
    
    
    
    hastighed = hastighedsgraf
    sted = stedgraf
    Nytid = acceleration[0]
    Nyacceleration = acceleration[1]
    middelacceleration = acceleration[1]
    
    MA = np.mean(middelacceleration)
    MAXH = np.max(hastighed)
    MINH = np.min(hastighed)
    STED = np.abs(np.abs(np.min(sted))-np.abs(np.max(sted)))
    
    print(f"Gennemsnitlig acceleration: {MA:.3f} m/s\u00B2 \n")
    print(f"Maksimal hastighed: {MAXH:.3f} m/s\n")
    print(f"Minimal hastighed: {MINH:.3f} m/s\n")
    print(f"Tilbagelagt afstand: {STED:.3f} m\n")
    
    ax[0].plot(acceleration[0],acceleration[1], linewidth = 1, color = "black", label = "Acceleration")
    
    
    ax[1].plot(Nytid,hastighed,linewidth = 1,color = "black", label = "Hastighed")
    
    ax[2].plot(Nytid,sted, linewidth = 1, color = "black", label = "Position")
    
    # ax[0].set_title(r"Acceleration (m/s$^2$)")
    ax[0].set_ylabel(r"Acceleration (m/s$^2$)") 
    
    # ax[1].set_title(r"Acceleration (m/s$^2$)")    
    # ax[2].set_title(r"Hastighed (m/s)")
    ax[1].set_ylabel(r"Hastighed (m/s)")
    
    # ax[3].set_title(r"sted (m)")
    ax[2].set_xlabel("tid (s)")
    ax[2].set_ylabel(r"sted (m)")
    for i in range(3):
        ax[i].legend()
    
    
    
    return 

def søg_bevægelsesligningerne(acceleration, fra, til):
    
    start = 0
    slut = 0
    
for i in range(len(acceleration[0])):
        if acceleration[0][0] >= fra:
            start = 0
            break
        if acceleration[0][i] >= fra:
            start += i
            break

    for i in range(len(acceleration[0])):
        if acceleration[0][-1] <= til:
            slut += -1
            break

        if acceleration[0][i] >= til:
            slut += i
            break
    
    hastighedsgraf = []
    stedgraf = []
    hastighed_global = 0
    sted_global = 0
    tid_forrige = acceleration[0][start]

    for i,j in zip(acceleration[1][start:slut],acceleration[0][start:slut]):
        dt = j-tid_forrige
        hastighed_global+=i*dt
        hastighedsgraf.append(hastighed_global)
        tid_forrige = j

    for i,j in zip(hastighedsgraf,acceleration[0][start:slut]):
        dt = j-tid_forrige
        sted_global+=i*dt
        stedgraf.append(sted_global)
        tid_forrige = j
    
    fig, ax = plt.subplots(4,1,figsize = (12,16))
    
    
    
    hastighed = hastighedsgraf
    sted = stedgraf
    Nytid = acceleration[0][start:slut]
    Nyacceleration = acceleration[1][start:slut]
    middelacceleration = acceleration[1][start+10:slut-10]
    
    MA = np.mean(middelacceleration)
    MAXH = np.max(hastighed)
    MINH = np.min(hastighed)
    STED = np.abs(np.abs(np.min(sted))-np.abs(np.max(sted)))
    
    print(f"Gennemsnitlig acceleration: {MA:.3f} m/s\u00B2 \n")
    print(f"Maksimal hastighed: {MAXH:.3f} m/s\n")
    print(f"Minimal hastighed: {MINH:.3f} m/s\n")
    print(f"Tilbagelagt afstand: {STED:.3f} m\n")
    
    ax[0].plot(acceleration[0],acceleration[1], linewidth = 1, color = "black")
    ax[0].plot([acceleration[0][start],acceleration[0][start]],[min(acceleration[1])+0.1*min(acceleration[1]),max(acceleration[1])+0.1*max(acceleration[1])],ls = '--',color = "#d62728",label = " Område der undersøges")
    ax[0].plot([acceleration[0][start],acceleration[0][slut]],[max(acceleration[1])+0.1*max(acceleration[1]),max(acceleration[1])+0.1*max(acceleration[1])],ls = '--',color = "#d62728")
    ax[0].plot([acceleration[0][slut],acceleration[0][slut]],[min(acceleration[1])+0.1*min(acceleration[1]),max(acceleration[1])+0.1*max(acceleration[1])],ls = '--',color = "#d62728")
    ax[0].plot([acceleration[0][start],acceleration[0][slut]],[min(acceleration[1])+0.1*min(acceleration[1]),min(acceleration[1])+0.1*min(acceleration[1])],ls = '--',color = "#d62728")

    
    ax[1].plot(Nytid,Nyacceleration, linewidth = 1, color = "black",label = "Acceleration" )
    ax[1].plot([acceleration[0][start+10],acceleration[0][slut-10]],[np.mean(middelacceleration),np.mean(middelacceleration)],ls = '--',color = "#d62728",label = f"Gennemsnit: {np.mean(middelacceleration):.5f} m/s\u00B2")
    ax[1].plot([Nytid[10],Nytid[10]],[np.min(middelacceleration),np.max(middelacceleration)],ls = '--',color = "#d62728")
    ax[1].plot([Nytid[-10],Nytid[-10]],[np.min(middelacceleration),np.max(middelacceleration)],ls = '--',color = "#d62728")
    
    
    ax[2].plot(Nytid,hastighed,linewidth = 1,color = "black", label = "Hastighed")
    
    ax[3].plot(Nytid,sted, linewidth = 1, color = "black", label = "Position")
    
    # ax[0].set_title(r"Acceleration (m/s$^2$)")
    ax[0].set_ylabel(r"Acceleration (m/s$^2$)") 
    
    # ax[1].set_title(r"Acceleration (m/s$^2$)")
    ax[1].set_ylabel(r"Acceleration (m/s$^2$)") 
    
    # ax[2].set_title(r"Hastighed (m/s)")
    ax[2].set_ylabel(r"Hastighed (m/s)")
    
    # ax[3].set_title(r"sted (m)")
    ax[3].set_xlabel("tid (s)")
    ax[3].set_ylabel(r"sted (m)")
    for i in range(4):
        ax[i].legend()
    
    
    
    return 

