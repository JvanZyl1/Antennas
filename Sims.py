'''
This is the program for the plotting of:
    a) The infinitesimal current element.
    b) The infinitesimal biconical antenna.
'''

import math
import numpy as np
import matplotlib.pyplot as plt

####Infinite current element####
rbylambda = np.arange(0.1, 100.1, 0.1)
C_c = 1 #This is "s"
Pow = []
for i in range(len(rbylambda)):
    Pc = -60 * (math.pi)**3 * C_c * math.cos(-4 * math.pi * rbylambda[i]) * rbylambda[i]
    Pow.append(Pc)

PowCE = False

if PowCE ==True:
    plt.figure(1)   
    plt.plot(rbylambda,Pow)
    plt.xlabel(r'$\frac{r}{\lambda}$')
    plt.ylabel("Power radiated [W] /s")
    plt.title("The power radiated for an infinitally long current element "
              "\n"
              r'$P = -60 \cdot \pi^3 \cdot s \cdot cos(-4 \pi \cdot \frac{r}{\lambda}) \cdot \frac{r}{\lambda}$')
    plt.show()

theta = np.arange(0.1, (math.pi + math.pi/20), (math.pi/10))
Zi = []
for i in range(len(theta)):
    inp = -120/((math.sin(theta[i]))**2)
    Zi.append(inp)

ZCE = False
if ZCE ==True:
    plt.figure(2)
    plt.plot(theta, Zi)
    plt.xticks(theta, ['0', 'π/10', 'π/5', '3π/10', '2π/5', 'π/2', '3π/5', '7π/10', '4π/5', '9π/10','π'])
    plt.margins(x=0)
    plt.ylabel("Input Impedance"
               r'$Z_i$')
    plt.xlabel(r'$\theta $ [rad]')
    plt.title("An infinitismal current element's input impedance"
              "\n"
              r'$Z_i = - \frac{120}{sin^2(\theta)}$')
    plt.show()



C_bi = 1 #H_0
# Using linspace so that the endpoint of 360 is included
actual = np.radians(np.linspace(0.1, 180.1, 20))
expected = np.arange(0.1, 100.1, 10)
 
lamb = 0.01 #m
dz = 1

r, theta = np.meshgrid(expected, actual)
values = np.ones((np.size(actual),np.size(expected)))
URe = np.ones((np.size(actual),np.size(expected)))
UIm = np.ones((np.size(actual),np.size(expected)))
Rr = np.ones((np.size(actual),np.size(expected)))
for i in range(len(actual)):
    for j in range(len(expected)):
        value =  120 * (math.pi**2) * (C_bi**2)*math.cos(4 * math.pi * expected[j-1])*(math.tan(actual[i-1]/2))**(-1)
        Re = (C_bi)**2 * np.log((math.tan(actual[i-1]))**(-1))*math.cos(-4*math.pi*expected[j-1])
        Im = (C_bi)**2 * np.log((math.tan(actual[i-1]))**(-1))*math.sin(-4*math.pi*expected[j-1])
        R_rad = - 60 * math.pi *(1/lamb)*math.sin(- expected[j-1]* math.pi * (1/lamb)) * math.sin(actual[i-1]) * dz**2
        Rr[i][j] = Rr[i][j]*R_rad
        URe[i][j] = URe[i][j]*Re
        UIm[i][j] = UIm[i][j]*Im
        values[i][j] = values[i][j]*value
        
plotIm = False
plotRe = False
plotP = False
plotR = True

if plotR == True:
    plt.figure(5)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    pc = ax.contourf(theta, r, Rr)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    plt.xlabel(r'$r$')
    plt.title("For an infinitesimal current element the radiation resistance."
              "\n"
              r'$R_{rad} = -\frac{60 \cdot \pi}{\lambda \cdot r} \cdot s^2 \cdot sin(-\beta \cdot r) \cdot sin(\theta)$'
              "\n"
              r'$\lambda = 1 mm, dz = 1$')
    cbar = plt.colorbar(pc)
    cbar.set_label(r'$R_{rad}$', rotation=270)
    plt.show()

if plotIm == True:
    plt.figure(5)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    pc = ax.contourf(theta, r, URe)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    plt.xlabel(r'$\frac{r}{\lambda}$')
    plt.title("For an ifinitesimal biconical antenna, the energy stored in the Real plane."
              "\n"
              r'$U = H_0^2 \cdot log|cot(\frac{\theta_{\frac{c}{2}}}{2})| \cdot [j \cdot sin(-4 \pi \cdot \frac{r}{\lambda})]$')
    cbar = plt.colorbar(pc)
    cbar.set_label(r'$U_{stored}$', rotation=270)
    plt.show()
if plotRe == True:
    plt.figure(5)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    pc = ax.contourf(theta, r, UIm)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    plt.xlabel(r'$\frac{r}{\lambda}$')
    plt.title("For an ifinitesimal biconical antenna, the energy stored in the Imaginary plane."
              "\n"
              r'$U = H_0^2 \cdot log|cot(\frac{\theta_{\frac{c}{2}}}{2})| \cdot [cos(-4 \pi \cdot \frac{r}{\lambda})]$')
    cbar = plt.colorbar(pc)
    cbar.set_label(r'$j \cdot U_{stored}$', rotation=270)
    plt.show()

if plotP == True:
    plt.figure(3)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    pc = ax.contourf(theta, r, values)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    plt.xlabel(r'$\frac{r}{\lambda}$')
    plt.title("For an ifinitesimal biconical antenna, the power radiated."
              "\n"
              r'$P = 120 \cdot \pi^2 \cdot H_0^2 \cdot cos(4 \cdot \pi \cdot \frac{r}{\lambda}) \cdot cot(\frac{\theta_{\frac{c}{2}}}{2})$')
    cbar = plt.colorbar(pc)
    cbar.set_label('Power Radiated [W]', rotation=270)
    plt.show()


theta = np.arange(0.1, (math.pi/2 + math.pi/20), (math.pi/20))
InputImpedance = []
for i in range(len(theta)):
    Zi = 120* math.log(abs(math.tan(theta[i]/2))**(-1))
    InputImpedance.append(Zi)
 
ZIBi = True
if ZIBi ==True:
    plt.figure(4)
    plt.plot(theta, InputImpedance)
    plt.xticks(theta, ['0', 'π/20', 'π/10', '3π/20', 'π/5', 'π/4', '3π/10', '7π/20', '2π/5', '9π/20', 'π/2'])
    plt.margins(x=0)
    plt.ylabel("Input Impedance"
               "\n"
               r'$Z_i$')
    plt.xlabel("Half-Cone Angle [rad]")
    plt.title("An infinitismal biconical antenna, the input impedance"
              "\n"
              r'$Z_i = 120 \cdot ln|cot(\frac{\theta_{\frac{c}{2}}}{2})|$')
    plt.show()