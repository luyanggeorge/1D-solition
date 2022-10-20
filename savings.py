from firedrake import *
import os.path
from pyop2.mpi import COMM_WORLD

"""
    *************************************************************
    *                        Saving Files                       *
    ************************************************************* """
def saving_files(save_path):
    save_waves = File(os.path.join(save_path, "waves.pvd"))
    save_WM = File(os.path.join(save_path, "Wavemaker.pvd"))
    Energy_file = open(os.path.join(save_path, 'energy.txt'), 'w')
    README_file = open(os.path.join(save_path, 'README.txt'), 'w')
    return save_waves, save_WM, Energy_file, README_file

""" *******************************************************
    *                  Compute the energy                 *
    ******************************************************* """

# yl update: XUVW. Delete: dWM_dt, WM, dWM_dy, x_coord
# add: b, C_ij, B_ij, FWF

def save_energy(dim, Lw, H0, g, H, h_n0, psi_1_n0, hat_psi_n0, Uu, Ww, VoW, IoW, WH, b, C11, CN1, CNN, B11, B1N, BN1, BNN, FWF,  
                A11, AN1, A1N, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN):
    if dim =="3D":
        energy = assemble(( 0.5*VoW*h_n0 * ((psi_1_n0.dx(0)**2)*M11 + dot(hat_psi_n0.dx(0), (2.0*MN1*psi_1_n0.dx(0)+dot(MNN,hat_psi_n0.dx(0)))))\
                           +0.5* Ww*h_n0 * ((psi_1_n0.dx(1)**2)*M11 + dot(hat_psi_n0.dx(1), (2.0*MN1*psi_1_n0.dx(1)+dot(MNN,hat_psi_n0.dx(1)))))\
                                +Uu*h_n0 * (psi_1_n0.dx(0)*(M11*psi_1_n0.dx(1) + dot(M1N,hat_psi_n0.dx(1))) \
                                                      + dot(hat_psi_n0.dx(0), (MN1*psi_1_n0.dx(1) +dot(MNN,hat_psi_n0.dx(1)))))\
                           -( VoW*h_n0.dx(0)+ Uu*h_n0.dx(1) ) * ( psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_n0)) \
                                                                 + dot(hat_psi_n0.dx(0), (DN1*psi_1_n0 + dot(DNN, hat_psi_n0))) ) \
                           -( Ww*h_n0.dx(1) + Uu*h_n0.dx(0) ) * ( psi_1_n0.dx(1)*(D11*psi_1_n0 + dot(D1N,hat_psi_n0)) \
                                                                 + dot(hat_psi_n0.dx(1), (DN1*psi_1_n0 + dot(DNN, hat_psi_n0)))) \
                           -H0*( VoW*b.dx(0)+ Uu*b.dx(1) ) * ( psi_1_n0.dx(0)*(B11*psi_1_n0 + dot(B1N,hat_psi_n0)) \
                                                                 + dot(hat_psi_n0.dx(0), (BN1*psi_1_n0 + dot(BNN, hat_psi_n0))) ) \
                           -H0*( Ww*b.dx(1) + Uu*b.dx(0) ) * ( psi_1_n0.dx(1)*(B11*psi_1_n0 + dot(B1N,hat_psi_n0)) \
                                                                 + dot(hat_psi_n0.dx(1), (BN1*psi_1_n0 + dot(BNN, hat_psi_n0)))) \
                           +     (1/(2.0*h_n0)) * (S11*psi_1_n0**2 + dot(hat_psi_n0, (2.0*SN1*psi_1_n0 + dot(SNN,hat_psi_n0))))*\
                                             ( VoW*h_n0.dx(0)**2 + Ww*h_n0.dx(1)**2 + 2*Uu*h_n0.dx(0)*h_n0.dx(1))
                           + (H0*H0/(2.0*h_n0)) * (A11*psi_1_n0**2 + dot(hat_psi_n0, (2.0*AN1*psi_1_n0 + dot(ANN,hat_psi_n0))))*\
                                             ( FWF*(VoW*b.dx(0)**2 + Ww*b.dx(1)**2 + 2.0*Uu*b.dx(0)*b.dx(1)) + Ww )\
                           +          (H0/h_n0) * (C11*psi_1_n0**2 + dot(hat_psi_n0, (2.0*CN1*psi_1_n0 + dot(CNN,hat_psi_n0))))*\
                                             ( VoW*b.dx(0)*h_n0.dx(0) + Uu*(b.dx(0)*h_n0.dx(1)+b.dx(1)*h_n0.dx(0)) + Ww*b.dx(1)*h_n0.dx(1))\
                           + 0.5*H0*g*Ww*h_n0*h_n0 - H0*g*WH*h_n0 )*dx) 

                           # yl notes: JCP paper Eq(17). full version!
    elif dim=="2D":
        energy = assemble(( 0.5*(h_n0*IoW)*(Lw**2)* ((psi_1_n0.dx(0)**2)*M11 \
                                                        + dot(hat_psi_n0.dx(0), (2.0*MN1*psi_1_n0.dx(0)+dot(MNN,hat_psi_n0.dx(0))))) \
                           -(IoW*(Lw**2)*h_n0.dx(0)) * (psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_n0))\
                                                        + dot(hat_psi_n0.dx(0),(DN1*psi_1_n0+dot(DNN, hat_psi_n0))))\
                           + (1/h_n0)*(0.5*IoW*(h_n0.dx(0)**2)*(Lw**2)) * (S11*psi_1_n0**2 + dot(hat_psi_n0, 2.0*SN1*psi_1_n0+dot(SNN,hat_psi_n0)))\
                           + 0.5*(Ww*H0*H0/h_n0) * (A11*psi_1_n0**2 + dot(hat_psi_n0, (2.0*AN1*psi_1_n0 + dot(ANN,hat_psi_n0))))\
                           + 0.5*H0*g*Ww*h_n0*h_n0 - H0*g*WH*h_n0 )*dx) # yl notes: JCP paper. MSA

                        # + 0.5*H0*g*Ww*(h_n0-H)**2)*dx) # FG's energy expression
    return energy


#----------------------------------------------------------------------#
#                        Surface solutions (x,y)                       #
#----------------------------------------------------------------------#
def x_to_xy(h_n0, WM, psi_1_n0, hat_psi_n0, h_2D, WM_2D, psi_s_2D, psi_i_2D, Indx):
    for i in range(len(h_n0.dat.data[:])):
        h_2D.dat.data[Indx[i]] = h_n0.dat.data[i]
        WM_2D.dat.data[Indx[i]] = WM.dat.data[i] # yl updated.
        psi_s_2D.dat.data[Indx[i]]=psi_1_n0.dat.data[i]
        psi_i_2D.dat.data[Indx[i],:] = hat_psi_n0.dat.data[i,:]
# yl notes: len(h_n0.dat.data[:]) = Nx+1
# h_2D.dat.data.shape = ((Nx+1)*2, )
# hat_psi_n0.dat.data.shape = (Nx+1,n_z), psi_i_2D.dat.data.shape = ((Nx+1)*2,n_z)
# extend to 2D: copy data from y=0 to y=Ly

#------------------------------------------------------------------------#
#                           3D solution (x,y,z)                          #
#------------------------------------------------------------------------#
def phi_projection(i, n_z, waves, Indz, psi_s, psi_i):
    if i==n_z:                                                   # if i=1,
        waves.dat.data[Indz[i]] = psi_s.dat.data[:]       # phi(z_i)=psi_1
    else:                                                        # if i>1,
        waves.dat.data[Indz[i]] = psi_i.dat.data[:,n_z-1-i] # phi(z_i)=psi_i
# yl notes: assembling phi, including the surface and the interior values
# yl question: how is data storaged in psi_i (based on a Vector Function Space)?

#-------------------------------------------------------------------------------------#
#                                 Transform the domain                                #
#-------------------------------------------------------------------------------------#

#____________________________________ z-transform ____________________________________#
def z_transform(mesh_3D, n_z, h_2D, beach_2D, H0, Indz):
    for i in range(0, n_z+1):                                          # for each layer
        mesh_3D.coordinates.dat.data[Indz[i],2]*=h_2D.dat.data[:]/H0   # z -> z*h/H0
        mesh_3D.coordinates.dat.data[Indz[i],2]+=beach_2D.dat.data[:]  # z -> z+b(x)

#____________________________________ x-transform ____________________________________#
def x_transform(mesh_3D, Lw, WM_3D, Indw):
    for i in range(0,len(Indw)): # x -> R + x*(Lw-R)/Lw
        mesh_3D.coordinates.dat.data[Indw[i],0]*=(Lw-WM_3D.dat.data[Indw[i]])/Lw
        mesh_3D.coordinates.dat.data[Indw[i],0]+=WM_3D.dat.data[Indw[i]]

def save_README(README_file, Lx, Ly, H0, xb, sb, res, Nx, Ny, n_z, gamma, Tw, w, t_stop, Lw, scheme, 
                dt, t, jours, heures, minutes, secondes, comp_time, before_it):
    README_file.write('______________________________________\n')
    README_file.write('                Summary               \n')
    README_file.write('--------------------------------------\n\n')

    README_file.write('------ Dimensions of the domain ------\n')
    README_file.write('Length Lx: %-5s m\n' %(str(Lx)))
    README_file.write('Length Ly: %-5s m\n' %(str(Ly)))
    README_file.write('Depth  H0: %-5s m\n' %(str(H0)))
    README_file.write('Beach start: %-5s m\n' %str(xb))
    README_file.write('Beach slope: %-5s\n\n' %str(sb))

    README_file.write('----------- Mesh resolution ----------\n')
    README_file.write('In x: %-5s m (%-5s elements)\n' % (str(res),str(Nx)))
    README_file.write('In y: %-5s m (%-5s elements)\n' % (str(res),str(Ny)))
    README_file.write('In z: %-5s m (%-5s elements)\n\n' % (str(H0/n_z),str(n_z)))

    README_file.write('-------------- Wavemaker -------------\n')
    README_file.write('Amplitude: %-5s m\n' %(str(gamma)))
    README_file.write('Period: %-5s s\n' %(str(Tw)))
    README_file.write('Frequency: %-5s /s\n' %(str(w)))
    README_file.write('Stops after %-5s periods ( = %-5s s)\n' %(str(t_stop/Tw), str(t_stop)))
    README_file.write('Lw: %-5s m\n\n' %(str(Lw)))

    README_file.write('--------------- Solver ---------------\n')
    if scheme=="SE":
        README_file.write('1st order Symplectic-Euler scheme\n\n')
    else:
        README_file.write('2nd order Stormer-Verlet scheme\n\n')

    README_file.write('------------- Final time -------------\n')
    README_file.write('Tend: %-5s s\n' %(str(t)))
    README_file.write('dt: %-5s s\n' %(str(dt)))
    README_file.write('Before time loop: %-5s s\n' %(str(before_it)))
    README_file.write('Computational time: %-5s d %-5s h %-5s m %-5s s\n' %(str(jours),str(heures),str(minutes), str(secondes)))
    README_file.write('That is in total: %-5s s\n' %(str(comp_time)))