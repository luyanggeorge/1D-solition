import numpy as np
import math
import matplotlib.pyplot as plt
import os.path

def get_data(file_path, file_name):
	file = os.path.join(file_path, file_name)
	xn = np.loadtxt(file, usecols=0)
	h  = np.loadtxt(file, usecols=1)
	phis = np.loadtxt(file, usecols=2)
	phib = np.loadtxt(file, usecols=3)
	return xn, h, phis, phib

def phi_expr(z):
    phi=a_phi*(Phi - 0.5*mu*np.power(z/H0,2)*Phi_2x + (mu*mu/24)*np.power(z/H0,4)*Phi_4x)
    return phi

save_path='data/t=T0+20.png'

g = 9.81
A = 0.1
H0= 20
e = 0.05
mu= e*e
a_eta=A*H0*e
a_phi=e*H0*np.sqrt(g*H0/mu)
c1 = 2*np.sqrt(4*e*A/3)

Lx = 8000
x0 = 0
x1 = Lx
dx = 1.0

speed= (1+0.5*A*e)*np.sqrt(g*H0)
dt = 0.02
T0 = round(0.5*Lx/speed)
Tend = T0 + 10

t = T0+10

path ='data/SE/2D/KPE_sech_2Aug/'
tt   = format(t,'.3f')
name = tt+'.txt'

xn, hn, phisn, phibn = get_data(path, name)

xi = (np.sqrt(3*e*A)/(2*H0)) * (xn - speed*t)
Phi = np.sqrt(4*e*A/3)*(np.tanh(xi)+1)
Phi_2x = -(A*e/mu) * np.sqrt(3*e*A) * np.tanh(xi) * np.power(np.cosh(xi),-2)
Phi_4x = 3 * np.power(A*e/mu,2) * np.sqrt(3*e*A) *\
                (2*np.tanh(xi)*np.power(np.cosh(xi),-4) - np.power(np.tanh(xi),3)*np.power(np.cosh(xi),-2)) 
    
eta = a_eta*np.power(np.cosh(xi),-2)
phi_b = phi_expr(0)
phi_s = phi_expr(H0+eta)

fig, (ax1, ax2, ax3) = plt.subplots(3,figsize=(9,9),constrained_layout=True)
ax1r= ax1.twinx()
ax2r= ax2.twinx() # compare the difference between numerical and exact solution
ax3r= ax3.twinx()

line1, = ax1.plot(xn, eta, 'k-',label='exact')
line1n,= ax1.plot(xn, hn-H0, 'r--',label='numerical')
line1d, = ax1r.plot(xn, eta-(hn-H0), 'm-.',label='diff')

line2, = ax2.plot(xn, phi_b, 'k-', label='exact')
line2n,= ax2.plot(xn, phibn, 'c--',label='numerical')
line2d, = ax2r.plot(xn, phi_b-phibn, 'm-.',label='diff')

line3, = ax3.plot(xn, phi_s, 'k-',label='exact',)
line3n,= ax3.plot(xn, phisn, 'c--',label='numerical')
line3d, = ax3r.plot(xn, phi_s-phisn, 'm-.',label='diff')

time_text=ax1.text( 0.8, 0.9, '', transform=ax1.transAxes, fontsize='large',
                   bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)) )

ax1.set_xlim(x0, x1)
#ax1.set_ylim(-0.1*a_eta, 1.1*a_eta)
ax1.set_title('$\eta (x,t)$',fontsize='x-large')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$\eta$')
ax1.grid()
ax1.legend(handles=[line1,line1n,line1d], loc='upper left', fontsize='large')

ax2.set_xlim(x0, x1)
ax2.set_ylim(-0.1*c1*a_phi, 1.1*c1*a_phi)
ax2.set_title('$\phi (x,z=0,t)$',fontsize='x-large')
ax2.set_xlabel('$x$')
ax2.set_ylabel('$\phi$')
ax2.grid()
ax2.legend(handles=[line2,line2n,line2d], loc='upper left', fontsize='large')

ax3.set_xlim(x0, x1)
ax3.set_ylim(-0.1*c1*a_phi, 1.1*c1*a_phi)
ax3.set_title('$\phi (x,z=h,t)$',fontsize='x-large')
ax3.set_xlabel('$x$')
ax3.set_ylabel('$\phi$')
ax3.grid()
ax3.legend(handles=[line3,line3n,line3d], loc='upper left', fontsize='large')

time_text.set_text('$t = %.2f$' %t)

plt.show()
#plt.savefig(save_path,dpi=300)