'''
This is the program for the plotting of:
    a) The finite current element.
'''

import math
import numpy as np
import matplotlib.pyplot as plt

freq = 100*(10**9)                  #Frequency of wave [Hz]
T = 1/freq                             #Time period [s]
I_0 = 2                             #Primary current applied [A]
l = 0.5                             #Length of cone along diagonal [m]   
epsilon0 = 8.8541878128 * (10**(-12)) # Vacuum Permeability [F⋅m−1]
omega = 2*math.pi*freq              #Angular velocity [rad/s] is 2pi/T = 2*pi*f
c = 3 * 10**8                       #Speed of light [m/s]

C_1 = ((I_0**2))/(12*math.pi*epsilon0)

t = np.linspace(0, 2*T+T/4, 4)
r = np.linspace(0.001, l-0.2, 50)
thetahc = (np.arcsin(r/l))
Pempty = np.ones((np.size(t),np.size(r),np.size(thetahc)))
X,Y = np.meshgrid(thetahc,r)
m = np.degrees(np.max(thetahc))

for i in range(len(t)):
    for j in range(len(thetahc)):
        for k in range(len(r)):
            Pow = C_1 * ((r[k]**2)/(math.sin(thetahc[j]))) * math.cos(4 * math.pi * freq * (t[i] - (r[k]/c))) * ((-4 * (math.pi**2)*(freq**2))/(c**2) + 1/(c*r[k]))
            Pempty[i,j,k] = Pempty[i,j,k]*Pow
plotpls = True
time = False
PT0 = Pempty[0,:,:]
PT1 = Pempty[1,:,:]
PT2 = Pempty[2,:,:]
PT3 = Pempty[3,:,:]

if plotpls == True:
    fig, ax = plt.subplots(2, 2, subplot_kw=dict(projection='polar'))
    CS = ax[0,0].contourf(X, Y, PT0)
    ax[0,0].set_title('t=0')
    ax[0,0].set_xlabel(r'R [m]')
    cbar = plt.colorbar(CS, ax = ax[0,0])
    cbar.set_label(r'$P[W]$', rotation=270)
    ax[0,0].set_thetamin(0)
    ax[0,0].set_thetamax(m)
    
    CSa = ax[1,0].contourf(X, Y, PT1)
    ax[1,0].set_title(r'$t = \frac{T}{4}$')
    ax[1,0].set_xlabel(r'R [m]')
    cbar = plt.colorbar(CSa, ax = ax[1,0])
    cbar.set_label(r'$P[W]$', rotation=270)
    ax[1,0].set_thetamin(0)
    ax[1,0].set_thetamax(m)
    
    CSb = ax[0,1].contourf(X, Y, PT2)
    ax[0,1].set_title(r'$t = \frac{T}{2}$')
    ax[0,1].set_xlabel(r'R [m]')
    cbar = plt.colorbar(CSb, ax = ax[0,1])
    cbar.set_label(r'$P[W]$', rotation=270)
    ax[0,1].set_thetamin(0)
    ax[0,1].set_thetamax(m)
    
    CSc = ax[1,1].contourf(X, Y, PT3)
    ax[1,1].set_title(r'$t = \frac{3T}{4}$')
    ax[1,1].set_xlabel(r'R [m]')
    cbar = plt.colorbar(CSc, ax = ax[1,1])
    cbar.set_label(r'$P[W]$', rotation=270)
    ax[1,1].set_thetamin(0)
    ax[1,1].set_thetamax(m)
a = r[25]
b = thetahc[25]
if time == True:
    plt.plot(t, Pempty[:,25,25])
    plt.xlabel("Time [s]")
    plt.ylabel("Power [W]")
    string = str("Power radiated [w], for a r = " + str(round(a, 3)) + r', $\theta$ = ' + str(round(b, 3)) + "and frequency = " + str(freq * 10**(-9)) + "GHz")
    plt.title(string)
    plt.show()