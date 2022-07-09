'''
This is the program for the plotting of:
    a) The finite biconical antenna
'''


import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
import math

gamma = 0.5772156649015328606065120 #Euler-Mascheroni constant
mu = 4*math.pi * 10**(-7) #Magnetic permeability of free space [H/m]

def saucy(angle,length):
    betalambda = 2 * math.pi
    l = length     #[wavelengths]
    betal = betalambda*l
    theta = angle
    Si2, Ci2 = sc.sici(2*betal)
    Si4, Ci4 = sc.sici(4*betal)
    Cin2 = gamma + np.log(2*betal) - Ci2
    
    Rm = 60*Cin2 + 30*(0.577 + np.log(betal) -2*Ci2 + Ci4)*np.cos(2*betal) + 30*(Si4 - 2*Si2)*np.sin(2*betal)
    Xm = 60*Si2 + 30*(Ci4 - np.log(betal) - 0.577)*np.sin(2*betal) - 30*(Si4)*np.cos(2*betal)
    Zk = 120* np.log(abs(np.tan(theta/2))**(-1))
    a = Zk - Xm*np.tan(betal)
    b = Rm * np.tan(betal)
    c = Rm
    d = Xm + Zk*np.tan(betal)
    ImZi = Zk*((a*c - b*d)/(c**2 + d**2))
    ReZi = Zk*((a*d + b*c)/(c**2 + d**2))
    ReL = mu* ((a*c - b*d) * (c**2 + d**2))/((a*c - b*d)**2 + (a*d + b*c)**2)
    ImL = mu* ((a*d - b*c) * (c**2 + d**2))/((a*c - b*d)**2 + (a*d + b*c)**2)
    return Rm, Xm, Zk, ImZi, ReZi, ReL, ImL

def length(y):
    betalambda = 2 * math.pi
    l = y     #[wavelengths]
    betal = betalambda*l
    Si2, Ci2 = sc.sici(2*betal)
    Si4, Ci4 = sc.sici(4*betal)
    Cin2 = gamma + np.log(2*betal) - Ci2
    Rm = 60*Cin2 + 30*(0.577 + np.log(betal) -2*Ci2 + Ci4)*np.cos(2*betal) + 30*(Si4 - 2*Si2)*np.sin(2*betal)
    Xm = 60*Si2 + 30*(Ci4 - np.log(betal) - 0.577)*np.sin(2*betal) - 30*(Si4)*np.cos(2*betal)
    return Rm, Xm

def currentsqr(radbylam, H_0):
    '''
    Outputs the square of the currents, and currents, in imaginary and real
    '''
    CI = 2 * math.pi * H_0
    CI2 = CI**2
    ReI = CI * np.cos(- 2*math.pi * radbylam)
    ImI = CI * np.sin(- 2*math.pi * radbylam)
    ReI2 = CI2 * 2 * np.cos(- 2*math.pi * radbylam) * np.sin(- 2*math.pi * radbylam)
    ImI2 = CI2 * np.cos(- 4*math.pi * radbylam)
    return ReI2, ImI2, ReI, ImI

def energy(ReL, ImL, ReI2, ImI2):
    ImU = 0.5 * (ReL * ReI2 - ImL * ImI2)
    ReU = 0.5 * (ReL * ImI2 - ReI2 * ImL)
    return ImU, ReU
    

x = np.radians(np.arange(-0, 90, 1))  #Theta [rad]
y = np.arange(0, 10.0, 0.001)   #l [lambda]
X, Y = np.meshgrid(x, y)

Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2

Rm, Xm, Zk, ImZi, ReZi, ReL, ImL = saucy(X,Y)
Rm1, Xm1 = length(y)

H_0 = 1

ReI2, ImI2, ReI, ImI = currentsqr(Y, H_0)
ImU, ReU = energy(ReL, ImL, ReI2, ImI2)


plotRm = True
plotImZi = True
plotEnergy = True

if plotRm == True:
    fig, ax = plt.subplots(1, 3)
    CSa = ax[0].plot(y, Rm1, 'b-')
    ax[0].set_xlabel("Length [Lambda]")
    ax[0].set_ylabel(r'$R_m$')
    CSb = ax[1].plot(y, Xm1, 'r-')
    ax[1].set_xlabel("Length [Lambda]")
    ax[1].set_ylabel(r'$X_m$')
    ax[2].plot(y, Rm1, 'b')
    ax[2].plot(y, Xm1, 'r-')
    ax[2].set_xlabel("Length [Lambda]")
    ax[2].set_ylabel(r'$Z_m$ component')
    ax[1].set_title(r'$Z_m = R_m + j \cdot X_m$')
    
if plotImZi == True:
    fig, ax = plt.subplots(1, 2)
    CS = ax[0].contourf(X, Y, ImZi)
    ax[0].set_title('Infinite Bi-Conical')
    ax[0].set_ylabel("Length [Lambda]")
    ax[0].set_xlabel(r'$\theta$')
    cbar = plt.colorbar(CS, ax = ax[0])
    cbar.set_label(r'$Im[Z_i]$', rotation=270)
    CS = ax[1].contourf(X, Y, ReZi)
    ax[1].set_title('Infinite Bi-Conical')
    ax[1].set_ylabel("Length [Lambda]")
    ax[1].set_xlabel(r'$\theta$')
    cbar = plt.colorbar(CS, ax = ax[1])
    cbar.set_label(r'$Re[Z_i]$', rotation=270)
    
if plotEnergy == True:
    fig, ax = plt.subplots(1, 2)
    CS = ax[0].contourf(X, Y, ImU)
    ax[0].set_title('Infinite Bi-Conical')
    ax[0].set_ylabel("Length [Lambda]")
    ax[0].set_xlabel(r'$\theta$')
    cbar = plt.colorbar(CS, ax = ax[0])
    cbar.set_label(r'$Im[U_{stored}]$', rotation=270)
    CS = ax[1].contourf(X, Y, ReU)
    ax[1].set_title('Infinite Bi-Conical')
    ax[1].set_ylabel("Length [Lambda]")
    ax[1].set_xlabel(r'$\theta$')
    cbar = plt.colorbar(CS, ax = ax[1])
    cbar.set_label(r'$Re[U_{stored}]$', rotation=270)