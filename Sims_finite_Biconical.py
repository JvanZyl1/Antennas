import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.special as sc

#https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sici.html



gamma = 0.5772156649015328606065120 #Euler-Mascheroni constant
actual = np.radians(np.linspace(0.1, 180.1, 20))
expected = np.arange(0.1, 100.1, 10)
r, theta = np.meshgrid(expected, actual)
Rm_l = np.ones((np.size(actual),np.size(expected)))
Xm_l = np.ones((np.size(actual),np.size(expected)))
Zk_l = np.ones((np.size(actual),np.size(expected)))
ImZi_l = np.ones((np.size(actual),np.size(expected)))
ReZi_l = np.ones((np.size(actual),np.size(expected)))
for i in range(len(actual)):
    for j in range(len(expected)):
        if i < 10:
            betalambda = 2 * math.pi
            l = actual[j-1]     #[wavelengths]
            betal = betalambda*l
            theta = expected[i-1]
            Si2, Ci2 = sc.sici(2*betal)
            Si4, Ci4 = sc.sici(4*betal)
            Cin2 = gamma + math.log(2*betal) - Ci2
            Rm = 60*Cin2 + 30*(0.577 + math.log(betal) -2*Ci2 + Ci4)*math.cos(2*betal) + 30*(Si4 - 2*Si2)*math.sin(2*betal)
            Xm = 60*Si2 + 30*(Ci4 - math.log(betal) - 0.577)*math.sin(2*betal) - 30*(Si4)*math.cos(2*betal)
            Zk = 120* math.log(abs(math.tan(theta/2))**(-1))
            a = Zk - Xm*math.tan(betal)
            b = Rm * math.tan(betal)
            c = Rm
            d = Xm + Zk*math.tan(betal)
            ImZi = Zk*((a*c - b*d)/(c**2 + d**2))
            ReZi = Zk*((a*d + b*c)/(c**2 + d**2))
            Rm_l[i,j] = Rm_l[i,j]*Rm
            Xm_l[i,j] = Xm_l[i,j]*Xm
            Zk_l[i,j] = Zk_l[i,j]*Zk
            ImZi_l[i,j] = ImZi_l[i,j]*ImZi
            ReZi_l[i,j] = ReZi_l[i,j]*ReZi

plotRm = True
plotXm = False
plotZk = False
plotImZi = False
plotReZi = False

print(np.shape(Rm_l))
print(np.shape(expected), np.shape(actual))

if plotRm == True:
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    pc = ax.contourf(theta, r, ImZi_l)
    plt.xlabel(r'$\frac{r}{\lambda}$')
    plt.title("For an ifinitesimal biconical antenna, the energy stored in the Real plane."
              "\n"
              r'$U = H_0^2 \cdot log|cot(\frac{\theta_{\frac{c}{2}}}{2})| \cdot [j \cdot sin(-4 \pi \cdot \frac{r}{\lambda})]$')
    cbar = plt.colorbar(pc)
    cbar.set_label(r'$U_{stored}$', rotation=270)
    plt.show()
    


