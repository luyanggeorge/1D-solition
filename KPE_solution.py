# single frame of the one-soliton exact solution of KPE
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

def phi_expr(z):
    phi=a_phi*(Phi - 0.5*mu*np.power(z/H0,2)*Phi_2x + (mu*mu/24)*np.power(z/H0,4)*Phi_4x)
    return phi

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
dx = 0.5

# To estimate the time step.
speed= (1+0.5*A*e)*np.sqrt(g*H0)
dt=dx/speed
print('dt_max = %.4f' %dt)

i=round(0.5*Lx/speed) # instant when x_peak=Lx/2
print(f't_start = {i}')

i=i+200 

fig, (ax1, ax2) = plt.subplots(2,figsize=(10,7),constrained_layout=True)
ax2r= ax2.twinx()
x = np.arange(x0, x1, dx)
line1, = ax1.plot([], [], 'b-')
line2, = ax2.plot([], [], 'k-',label='$z=0$')
line3, = ax2.plot([], [], 'c--',label='$z=H_0+\eta$')
line4, = ax2r.plot([], [], 'r--',label='diff')
time_text=ax1.text( 0.8, 0.9, '', transform=ax1.transAxes, fontsize='large',
                   bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)) )

ax1.set_xlim(x0, x1)
ax1.set_ylim(-0.1*a_eta, 1.1*a_eta)
ax1.set_title('$\eta (x,t)$',fontsize='x-large')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$\eta$')
ax1.grid()

ax2.set_xlim(x0, x1)
ax2.set_ylim(-0.1*c1*a_phi, 1.1*c1*a_phi)
ax2.set_title('$\phi (x,z,t)$',fontsize='x-large')
ax2.set_xlabel('$x$')
ax2.set_ylabel('$\phi$')
ax2.grid()
ax2.legend(handles=[line2,line3,line4], loc='upper left', fontsize='large')

xi = (np.sqrt(3*e*A)/(2*H0)) * (x - (1+0.5*A*e)*np.sqrt(g*H0)*i)
Phi = np.sqrt(4*e*A/3)*(np.tanh(xi)+1)
Phi_2x = -(A*e/mu) * np.sqrt(3*e*A) * np.tanh(xi) * np.power(np.cosh(xi),-2)
Phi_4x = 3 * np.power(A*e/mu,2) * np.sqrt(3*e*A) *\
                (2*np.tanh(xi)*np.power(np.cosh(xi),-4) - np.power(np.tanh(xi),3)*np.power(np.cosh(xi),-2)) 
    
eta = a_eta*np.power(np.cosh(xi),-2)
phi_b = phi_expr(0)
phi_s = phi_expr(H0+eta)
line1.set_data(x, eta)
line2.set_data(x, phi_b)
line3.set_data(x, phi_s)
line4.set_data(x, phi_b-phi_s)
time_text.set_text('$t = %.1f$' %i)

#print(np.amax(phi_b-phi_s))
plt.show()