import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

scheme='SE'
#save_path='data/h0025_TC2 vs TC2b.png'

file_1='data/'+scheme+'/2D/KPE_sech/energy.txt'


time1   =np.loadtxt(file_1, usecols=0)
en_tot1 =np.loadtxt(file_1, usecols=1)
h_1     =np.loadtxt(file_1, usecols=2)
psi_s_1  =np.loadtxt(file_1, usecols=3)
psi_b_1    =np.loadtxt(file_1, usecols=4)

en1 = []
for i in range(len(en_tot1)):
    en1.append(en_tot1[i]-en_tot1[0])

fig, [ax2, ax3]=plt.subplots(2)
fig.set_size_inches(9,9)
fig.set_tight_layout(True)

ax2.set_title('Energy variations', fontsize=20)
ax2.set_xlabel('Time [s]',fontsize=14)
ax2.set_ylabel('E(t) [J]',fontsize=14)
ax2.plot(time1,en1,'b-')
ax2.grid()

ax3.set_title('Water depth at (0,0)', fontsize=20)
ax3.set_xlabel('Time [s]',fontsize=14)
ax3.set_ylabel('h(x=0,y=0,t) [m]',fontsize=14)
ax3.plot(time1,h_1,'b-')
ax3.grid()

plt.show()
#plt.savefig(save_path,dpi=300)