import pdb
import time
import numpy as np
import os.path
from scipy.interpolate import interp1d
from vertical_discr_full import *
from savings import *
from firedrake import *
from firedrake.petsc import PETSc
import solvers_full as YL_solvers
from settings import *

start_time = time.perf_counter()

"""
    ****************************************
    *               Settings               *
    **************************************** """

input_data, scheme, dim, save_path, bottom, FWF, save_pvd = test_case()
H0, xb, sb, H_expr, Hend, Lx, Ly, Lw, res_x, res_y, n_z = domain(bottom)
g, gamma, speed, WM_expr, dWM_dt_expr, dWM_dy_expr, xi_expr, h_ex_expr, Phi_expr,\
Phi_2x_expr, Phi_4x_expr, phi_ex_expr = wavemaker(dim, H0, Ly, Lw, input_data)
T0, t, dt, Tend, t_stop, dt_save = set_time(speed, Lx, res_x, Ly, res_y, H0, n_z)

PETSc.Sys.Print('...settings loaded!')

PETSc.Sys.Print('Creation of the mesh across %d processes...' % COMM_WORLD.size)
    
"""
    ****************************************************
    *               Definition of the mesh             *
    **************************************************** """

#_________________ Vertical discretization ________________#
Nz = n_z+1         # Number of point in one vertical element

#________________ Horizontal discretization _______________#
Nx = round(Lx/res_x)    # Number of elements in x (round to the nearest integer)
Ny = round(Ly/res_y)    # Number of elements in y

#___________________________ Mesh _________________________#
if dim=="2D":                                   #(x,z)-waves
    # Generate a uniform mesh of an interval (L_x).
    hor_mesh = IntervalMesh(Nx,Lx) 
else:                                         #(x,y,z)-waves
    # Generate a rectangular mesh.
    hor_mesh = RectangleMesh(Nx,Ny,Lx,Ly,quadrilateral=True)
    # quadrilateral – (optional), creates quadrilateral mesh, defaults to False

PETSc.Sys.Print('...mesh created!')

PETSc.Sys.Print('Definition of the function...')

"""
    *************************************************
    *       Definition of the function spaces       *
    ************************************************* """
#___________________ For h and psi_1 ___________________#
V = FunctionSpace(hor_mesh, "CG", 1)
#_____________________ For hat_psi _____________________#
Vec = VectorFunctionSpace(hor_mesh, "CG", 1,dim=n_z)
# We might want the number of components in the vector to differ from the geometric dimension of the mesh. 
# We can do this by passing a value for the dim argument to the VectorFunctionSpace() constructor.
#_________________ Mixed function space ________________#
V_mixed = V*Vec # to solve simultaneous weak formulations

"""
    ******************************************************
    *            Definition of the functions             *
    ****************************************************** """

if scheme=="SE": #_________ Symplectic-Euler scheme _________#
    #______________________ At time t^n _____________________#
    h_n0 = Function(V)                                   # h^n
    psi_1_n0 = Function(V)                           # psi_1^n
    hat_psi_n0 = Function(Vec)                     # hat_psi^n

    #________________ At time t^{n+1} and t^* _______________#
    psi_1_n1 = Function(V)                       # psi_1^{n+1}
    w_n1 = Function(V_mixed)
    h_n1, hat_psi_star = split(w_n1)      # h^{n+1}, hat_psi^*
    hat_psi_n1 = Function(Vec)    # to visualise hat_psi^{n+1}
else: #________________ Stormer-Verlet scheme _______________#
    #______________________ At time t^n _____________________#
    h_n0 = Function(V)                                   # h^n
    psi_1_n0 = Function(V)                           # psi_1^n
    hat_psi_n0 = Function(Vec)                     # hat_psi^n

    #_______________ At time t^{n+1/2} and t^* ______________#
    w_half = Function(V_mixed)        # to obtain psi^{n+1/2},
    psi_1_half, hat_psi_star = split(w_half)   # and hat_psi^*

    #_______________ At time t^{n+1} and t^** _______________#
    psi_1_n1 = Function(V)                       # psi_1^{n+1}
    w_n1 = Function(V_mixed)              # to obtain h^{n+1},
    h_n1, hat_psi_aux = split(w_n1)         # and hat_psi^{**}
    hat_psi_n1 = Function(Vec)    # to visualise hat_psi^{n+1}


#_______________________ x coordinate _______________________#
# yl update
x = SpatialCoordinate(hor_mesh)
x_coord = Function(V).interpolate(x[0])

#________________ Beach topography b(x,y)____________________#
b = Function(V)                                       # b(x,y)

#_______________________ Depth at rest ______________________#
H = Function(V)                                         # H(x)

#_________________________ Wavemaker ________________________#
WM = Function(V)                                  # R(x,y;t^n)
dWM_dt = Function(V)                               # (dR/dt)^n
dWM_dy = Function(V)                               # (dR/dy)^n
if scheme=="SV":                         # For Stormer-Verlet:
    WM_half = Function(V)                   # R(x,y;t^{n+1/2})
    dWM_half_dt = Function(V)                # (dR/dt)^{n+1/2}
    dWM_half_dy = Function(V)                # (dR/dy)^{n+1/2}
WM_n1 = Function(V)                           # R(x,y;t^{n+1})
dWM_n1_dt = Function(V)                        # (dR/dt)^{n+1}
dWM_n1_dy = Function(V)                        # (dR/dy)^{n+1}

#______________________ Trial functions _____________________#
psi_1 = TrialFunction(V)      # psi_1^{n+1} for linear solvers
hat_psi = TrialFunction(Vec)# hat_psi^{n+1} for linear solvers

#_______________________ Test functions _____________________#
delta_h = TestFunction(V)                         # from dH/dh
delta_hat_psi = TestFunction(Vec)           # from dH/dhat_psi
w_t = TestFunction(V_mixed)                # from dH/dpsi_1...
delta_psi, delta_hat_star = split(w_t)    # ...and dH/dhat_psi
# yl update:
if scheme=="SV": 
    w_t_sv = TestFunction(V_mixed)
    delta_h_sv, delta_hat_psi_sv = split(w_t_sv)

phii_z = Function(V) # for initialising standing wave solution

xi  = Function(V)
Phi = Function(V)
Phi_2x = Function(V)
Phi_4x = Function(V)
z_bt =Function(V)

# ------SE------
# step1: use delta_psi and delta_hat_star to solve simultaneously for h^n+1 and psi_hat^*
# step2: use delta_h to solve for psi_1^(n+1)
# step3: use delta_hat_psi to update psi_hat^(n+1) using Laplace Eq.
# ------SV------
# step1: use delta_h_sv and delta_hat_psi_sv to solve simultaneously for psi_1^half and psi_hat^*
# step2: use delta_psi and delta_hat_star to solve simultaneously for h^n+1 and psi_hat^**
# step3: use delta_h to solve for psi_1^(n+1)
# step4: use delta_hat_psi to update psi_hat^(n+1) using Laplace Eq.

# yl update: XUVW
Xx=Function(V) 

Uu=Function(V)
Ww=Function(V)
Vv=Function(V)
VoW=Function(V)
IoW=Function(V)
WH=Function(V)  # W x H
XRt=Function(V) # X x dR/dt

Ww_n1=Function(V) # for SE and SV

if scheme=="SV":  # for SV
    Uu_half =Function(V)
    Ww_half =Function(V)
    Vv_half =Function(V)
    VoW_half=Function(V)
    IoW_half=Function(V)
    WH_half =Function(V)  # W x H
    XRt_half=Function(V)  # X x dR/dt

PETSc.Sys.Print('...functions created!')

PETSc.Sys.Print('Initalisation of the functions...')
"""
    ***********************************************************************************
    *                          Initialisation of the Functions                        *
    ***********************************************************************************"""
#---------------------------- Topography ----------------------------#             
# yl update
H.assign(H_expr(V,x))                             # Depth at rest H(x)
b.assign(H0-H)                                    # Beach b(x,y)

#----------------------------------------------------------------------------------------#
#                                       Wavemaker                                        #
#----------------------------------------------------------------------------------------#
if input_data=="measurements":#------------- Interpolate measurements -------------#
    # yl added
    t_data, wm_data, wm_vel_data = load_wavemaker(dt)
    wm_inter     = interp1d(t_data,wm_data)
    wm_vel_inter = interp1d(t_data,wm_vel_data)
    # yl updated
    WM_expr(WM, x, t, wm_inter)
    dWM_dt_expr(dWM_dt, x, t, wm_vel_inter)
    dWM_dy_expr(dWM_dy)
else:
    # yl update
    WM_expr(WM,x,t)                  # \tilde{R}(x,y;t)
    dWM_dt_expr(dWM_dt,x,t)          # d\tilde{R}/dt  
    dWM_dy_expr(dWM_dy,x,t)          # d\tilde{R}/dy

#----------------------------------------------------------------------------------------#
#                               Solutions Initialization                                 #
#----------------------------------------------------------------------------------------#
#____________________________ Initialization of Depth ___________________________________#
xi_expr(xi,x,t)
h_ex_expr(h_n0,xi)
w_n1.sub(0).assign(h_n0)

#h_n0.assign(H)        # h(x,y;t=0) = H(x)                                                  
#w_n1.sub(0).assign(H) # Extract the ith sub Function of this Function. In this case, h^{n+1}.

#_____________________ Velocity pot. at the surface: phi(x,y,z=h;t) _____________________#
Phi_expr(Phi,xi)
Phi_2x_expr(Phi_2x,xi)
Phi_4x_expr(Phi_4x,xi)

phi_ex_expr(psi_1_n0, Phi, Phi_2x, Phi_4x, h_n0)

#psi_1_n0.assign(0.0)                                         # \psi_1(x,y;t=0) = 0

#_____________________ Velocity pot. in depth: phi(x,y,z<h;t) _____________________#
for i in range(0,n_z):
    z=(H0/n_z)*(n_z-i-1) # after coordinate transfermation hat_z
    z_bt.assign((z/H0)*h_n0) # before coordinate transfermation
    phi_ex_expr(phii_z, Phi, Phi_2x, Phi_4x, z_bt)
    hat_psi_n0.dat.data[:,i] = phii_z.dat.data   # psi_i^n
    w_n1.sub(1).dat.data[:,i] = phii_z.dat.data  # psi_i^{*}

#for i in range(0,n_z):
#    hat_psi_n0.dat.data[:,i] = 0.0 # psi_i^n
#    w_n1.sub(1).dat.data[:,i] = 0.0  # psi_i^{*}

#______________________________ Stormer-Verlet scheme _____________________________#
# initialise psi_1_half and hat_psi_star for SV

if scheme=="SV": 
    w_half.sub(0).assign(psi_1_n0)
    for i in range(0,n_z):
        z=(H0/n_z)*(n_z-i-1)
        z_bt.assign((z/H0)*h_n0)
        phi_ex_expr(phii_z, Phi, Phi_2x, Phi_4x, z_bt)
        w_half.sub(1).dat.data[:,i]= phii_z.dat.data

#if scheme=="SV": 
#    w_half.sub(0).assign(0.0)
#    for i in range(0,n_z):
#        w_half.sub(1).dat.data[:,i]= 0.0

# yl update: XUVW
Xx.assign(x_coord-Lw)
Uu.assign(Xx*dWM_dy)
Ww.assign(Lw-WM)
Vv.assign(Lw*Lw+Uu*Uu)
VoW.assign(Vv/Ww)
IoW.assign(1/Ww)
WH.assign(Ww*H)
XRt.assign(Xx*dWM_dt)

if scheme=='SV' and input_data=='measurements':
    Uu_half.assign(0.0)
    Vv_half.assign(Lw*Lw)

PETSc.Sys.Print('...functions initialised!')

PETSc.Sys.Print('Assembling z-matrices...')

"""
    ************************
    * Compute the matrices *
    ************************ """
#_______ Initialization ______#
A = np.zeros((Nz,Nz))
B = np.zeros((Nz,Nz)) # FWF
C = np.zeros((Nz,Nz)) # FWF
M = np.zeros((Nz,Nz))
D = np.zeros((Nz,Nz))
S = np.zeros((Nz,Nz))
Ik = np.zeros((Nz,1))

#____ Filling the matrices ___#
for i in range(0,Nz):
    for j in range(0,Nz):
        A[i,j]=A_ij(i,j,n_z,H0)
        M[i,j]=M_ij(i,j,n_z,H0)
        D[i,j]=D_ij(i,j,n_z,H0)
        S[i,j]=S_ij(i,j,n_z,H0)
    Ik[i] = I_i(i,n_z,H0)

#________ Submatrices ________#
A11 = A[0,0]
A1N = as_tensor(A[0,1:])
AN1 = as_tensor(A[1:,0])
ANN = as_tensor(A[1:,1:])

M11 = M[0,0]
M1N = as_tensor(M[0,1:])
MN1 = as_tensor(M[1:,0])
MNN = as_tensor(M[1:,1:])

# yl comments:
# M1N or MN1: class 'ufl.tensors.ListTensor', len(M1N)=n_z, shape(M1N)=shape(MN1)=(n_z,), M1N[i]=MN1[i]
# MNN: class 'ufl.tensors.ListTensor', shape(MNN)=(n_z,n_z)

D11 = D[0,0]
D1N = as_tensor(D[0,1:])
DN1 = as_tensor(D[1:,0])
DNN = as_tensor(D[1:,1:])
# D1N[i]!=DN1[i]

S11 = S[0,0]
S1N = as_tensor(S[0,1:])
SN1 = as_tensor(S[1:,0])
SNN = as_tensor(S[1:,1:])

I1 = Ik[0,0]
IN=as_tensor(Ik[1:,0])

# yl added: full weak forms
if FWF==1:
    for i in range(0,Nz):
        for j in range(0,Nz):
            B[i,j]=B_ij(i,j,n_z,H0)
            C[i,j]=C_ij(i,j,n_z,H0)

B11 = B[0,0]
B1N = as_tensor(B[0,1:])
BN1 = as_tensor(B[1:,0])
BNN = as_tensor(B[1:,1:])

C11 = C[0,0]
C1N = as_tensor(C[0,1:])
CN1 = as_tensor(C[1:,0])
CNN = as_tensor(C[1:,1:])

PETSc.Sys.Print('... end of assembling!')

PETSc.Sys.Print('Initialisation of the solvers...')

"""
    ************************************************************************************************************************
    *                                                   Weak Formulations                                                  *
    ************************************************************************************************************************ """

if scheme=="SE": #_____________________________________________ Symplectic-Euler ______________________________________________#
    # yl update: XUVW, full 3D WF
    #------------------------ Step 1 : Update h at time t^{n+1} and psi_i at time t^* simulataneously: ------------------------#
    WF_h_psi = YL_solvers.WF_h_SE(dim, n_z, g, H0, Lw, dWM_dt, dt, delta_psi, delta_hat_star, h_n0, h_n1, psi_1_n0, hat_psi_star,
                                  Uu, Ww, VoW, IoW, WH, XRt, b, C11, CN1, CNN, B11, B1N, BN1, BNN, FWF,
                                  M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)
    # full 3D WF.
    #----------------------------------------- Step 2 : Update psi_1 at time t^{n+1}: -----------------------------------------#
    A_psi_s, L_psi_s = YL_solvers.WF_psi_SE(dim, g, H0, Lw, dWM_dt, dt, delta_h, psi_1, psi_1_n0, hat_psi_star, h_n1, 
                                            Uu, Ww, Ww_n1, VoW, IoW, WH, XRt, b, C11, CN1, CNN, B11, B1N, BN1, BNN, FWF,
                                            M11, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)

    #----------------------------------------- Step 3 : Update psi_i at time t^{n+1}: -----------------------------------------#
    A_hat, L_hat = YL_solvers.WF_hat_psi_SE(dim, H0, n_z, Lw, dWM_dt, dt, delta_hat_psi, hat_psi, h_n0, psi_1_n0, 
                                            Uu, Ww, VoW, IoW, M11, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)

elif scheme=="SV":#______________________________________________ Stormer-Verlet ______________________________________________#
    # yl update: XUVW, full 3D WF
    #--------------------------------------- Step 1 : Update psi_1^{n+1/2} and psi_i^*: ---------------------------------------#
    WF_psi_star = YL_solvers.WF_psi_half_SV(dim, n_z, g, H0, Lw, dt, delta_h_sv, delta_hat_psi_sv, psi_1_n0, psi_1_half, hat_psi_star, h_n0, dWM_half_dt, 
                                            Ww, Uu_half, VoW_half, Ww_half, XRt_half, WH_half, IoW_half, b, C11, CN1, CNN, B11, B1N, BN1, BNN, FWF,
                                            M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)
    # full 3D WF
    #----------------------------- Step 2 : Update h^{n+1} and psi_i at time t^** simulataneously: ----------------------------#
    WF_h_psi = YL_solvers.WF_h_SV(dim, n_z, Lw, H0, g, dt, delta_psi, delta_hat_star, h_n0, h_n1, psi_1_half, hat_psi_star, hat_psi_aux, 
                                  dWM_half_dt, Ww_half, Uu_half, VoW_half, XRt_half, IoW_half, b, C11, CN1, CNN, B11, B1N, BN1, BNN, FWF,
                                  M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)
    # full 3D WF
    #----------------------------------------- Step 3 : Update psi_1 at time t^{n+1}: -----------------------------------------#
    a_psi_1, L_psi_1 = YL_solvers.WF_psi_n1_SV(dim, H0, g, delta_h, Lw, dt, psi_1_half, psi_1, dWM_half_dt, hat_psi_aux, h_n1, 
                                               Ww_n1, Ww_half, Uu_half, VoW_half, XRt_half, WH_half, IoW_half, b, C11, CN1, CNN, B11, B1N, BN1, BNN, FWF,
                                               M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)

    #----------------------------------------- Step 4 : Update psi_i at time t^{n+1}: -----------------------------------------#
    A_hat, L_hat = YL_solvers.WF_hat_psi_SV(dim, n_z, Lw, H0, dt, delta_hat_psi, hat_psi, h_n0, psi_1_n0, 
                                            dWM_dt, Uu, Ww, VoW, IoW,
                                            M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)
"""
    **************************************************************************************
    *                                 Define the solvers                                 *
    ************************************************************************************** """

#____________________________________ Solvers parameters ____________________________________#
'''
param_h       = {'ksp_converged_reason':None, 'pc_type': 'fieldsplit','pc_fieldsplit_type': 'schur','pc_fieldsplit_schur_fact_type': 'upper'}               
param_psi     = {'ksp_converged_reason':None, 'ksp_type': 'preonly', 'pc_type': 'lu'}
param_hat_psi = {'ksp_converged_reason':None, 'ksp_type': 'preonly', 'pc_type': 'lu'}
'''
# before optimisation
param_h       = {'ksp_converged_reason':None}               
param_psi     = {'ksp_converged_reason':None}
param_hat_psi = {'ksp_converged_reason':None}

# By default, the solve call will use GMRES using an incomplete LU factorisation to precondition the problem.
# We may solve the system directly by computing an LU factorisation of the problem. 
# To do this, we set the pc_type to 'lu' and tell PETSc to use a “preconditioner only” Krylov method.
# see https://www.firedrakeproject.org/solving-interface.html#solving-linear-systems
# 'ksp_view': None => print the solver parameters
# 'snes_view': None => print the residual
# 'snes_monitor': None

#--------------------------------------------------------------------------------------------#
#                                      Symplectic-Euler                                      #
#____________________________________________________________________________________________#
if scheme=="SE":
    #_______________________ Variational solver for h (and hat_psi^*) _______________________#
    h_problem = NonlinearVariationalProblem(WF_h_psi, w_n1)
    h_solver = NonlinearVariationalSolver(h_problem, options_prefix="h_dt_imp", solver_parameters=param_h)

    #_____________________________ Variational solver for psi_1 _____________________________#
    psi_problem = LinearVariationalProblem(A_psi_s, L_psi_s, psi_1_n1)
    # yl comment:
    # In this linear solver the trial function is psi_1.
    # psi_1_n1 is a function holding the solution, or we place the solution in psi_1_n1.
    psi_solver = LinearVariationalSolver(psi_problem, options_prefix="psi1_dt_exp", solver_parameters=param_psi)

    #____________________________ Variational solver for hat_psi ____________________________#
    hat_psi_problem = LinearVariationalProblem(A_hat, L_hat, hat_psi_n0)
    # yl comment:
    # psi_1_n1 was created and initialised (for calculating energy at t=0), in this linear solver the trial function is hat_psi.
    hat_psi_solver = LinearVariationalSolver(hat_psi_problem, options_prefix="hat_psi_exp", solver_parameters=param_hat_psi)

#--------------------------------------------------------------------------------------------#
#                                       Stormer-Verlet                                       #
#____________________________________________________________________________________________#
if scheme=="SV":
    #_______________________ Variational solver for psi_1^{n+1/2} (and hat_psi^*) _______________________#
    psi_half_problem = NonlinearVariationalProblem(WF_psi_star, w_half)
    psi_half_solver = NonlinearVariationalSolver(psi_half_problem, options_prefix="psi1_dt2_imp", solver_parameters=param_h)
    
    #____________________________ Variational solver for h^{n+1} (and hat_psi^**) _______________________#
    h_problem = NonlinearVariationalProblem(WF_h_psi, w_n1)
    h_solver = NonlinearVariationalSolver(h_problem, options_prefix="h_dt_imp", solver_parameters=param_h)
    
    #_______________________ Variational solver for psi_1^{n+1} _______________________#
    psi_n1_problem = LinearVariationalProblem(a_psi_1, L_psi_1, psi_1_n1)
    psi_n1_solver = LinearVariationalSolver(psi_n1_problem, options_prefix="psi1_dt_exp", solver_parameters=param_psi)
    
    #____________________________ Variational solver for hat_psi ____________________________#
    hat_psi_problem = LinearVariationalProblem(A_hat, L_hat, hat_psi_n0)
    hat_psi_solver = LinearVariationalSolver(hat_psi_problem, options_prefix="hat_psi_exp", solver_parameters=param_hat_psi)

PETSc.Sys.Print('...solvers initialised!')

"""
    *************************************************************
    *                        Saving Files                       *
    ************************************************************* """
save_waves, save_WM, Energy_file, README_file = saving_files(save_path)

"""
    ****************************************************************************
    *                                 Saving mesh                              *
    ****************************************************************************"""
#---------------------------------------------------------------------------------#
#                      Save waves in the 3D free-surface domain                   #
#---------------------------------------------------------------------------------#

if dim=='2D': # Extend the 1D horizontal mesh (x) to 2D horizontal mesh (x,y)
    mesh_2D = RectangleMesh(Nx,1,Lx,Ly,quadrilateral=True)        # 2D surface mesh
    V_2D = FunctionSpace(mesh_2D,"CG",1)                  # 2D surface funct. space
    Vec_2D = VectorFunctionSpace(mesh_2D,"CG",1, dim=n_z)  # 2D vector funct. space
    h_2D = Function(V_2D)                                                  # h(x,y)
    psi_s_2D = Function(V_2D)                                         # psi_1 (x,y)
    psi_i_2D = Function(Vec_2D)                                       # psi_i (x,y)
    # yl updated:
    WM_2D = Function(V_2D)
    x2 = SpatialCoordinate(mesh_2D)                                       
    beach_s_2D = Function(V_2D).interpolate(conditional(le(x2[0],xb),0.0,sb*(x2[0]-xb)))
    # Extend the surface mesh in depth to obtain {0<x<Lx; 0<y<Ly; 0<z<H0}
    mesh_3D = ExtrudedMesh(mesh_2D,                   # horizontal mesh to extrude;
                           n_z,               # number of elements in the vertical;
                           layer_height=H0/(n_z),         # length of each element;
                           extrusion_type='uniform')     # type of extruded coord.;

else:# If the solutions are already (x,y)-dependent, we extend the domain in depth:
    mesh_3D = ExtrudedMesh(hor_mesh,                  # horizontal mesh to extrude;
                           n_z,               # number of elements in the vertical;
                           layer_height=H0/(n_z),         # length of each element;
                           extrusion_type='uniform')     # type of extruded coord.;

"""
    *****************************
    *      Function to save     *
    ***************************** """
#__________ Function Space _________#
V_3D = FunctionSpace(mesh_3D, "CG",1)
#____________ Functions ____________#
waves = Function(V_3D,name= "phi")  # yl notes: store all the phi(x,y,z,t^n0) in the domain.
WM_3D = Function(V_3D,name = "WM") 


"""
    **************************************************************************
    *                         Mapping and transforms                         *
    **************************************************************************"""
if dim=="2D":
    # Indices to map h(x) and phi(x) to h(x,y) and phi(x,y) :
    Indx = []
    nodes=len(hor_mesh.coordinates.dat.data)
    if nodes<=1001:
        for j in range(nodes):
            Indx.append([y for y in range(len(mesh_2D.coordinates.dat.data[:,0]))\
            if mesh_2D.coordinates.dat.data[y,0]==hor_mesh.coordinates.dat.data[j]])
    else:
        for j in range(nodes-2):
            Indx.append([2*(nodes-1-j),2*(nodes-1-j)+1])
        Indx.append([0,1])
        Indx.append([2,3])
# yl notes: len(hor_mesh.coordinates.dat.data) = Nx+1, len(Indx) = Nx+1, each element is a pair of numbers.
# For the jth set Indx[j], each element represents the No. of the node on the vertical line x=(j-1)*(Lx/Nx),j=1,2,...,Nx+1

# Index used to differentiate each vertical layer
Indz = []
for i in range(0,n_z+1):
    Indz.append([zz for zz in range(len(mesh_3D.coordinates.dat.data[:,2])) \
     if mesh_3D.coordinates.dat.data[zz,2] == mesh_3D.coordinates.dat.data[i,2]])
# yl notes: len(Indz)=n_z+1, each element is also a list (set) with (Nx+1)*(Ny+1) numbers.
# For the ith set Indz[i], each element represents the No. of the node on the plane z=(i-1)*(H0/n_z),i=1,2,...,n_z+1

# Index of the 3D funct. for which x<Lw. This is used to transform the 3D domain
# in x, to get back to the moving domain:
Test_x_Lw=Function(V_3D)
# yl updates:
x1 = SpatialCoordinate(mesh_3D) # Need to define the coordinates again. yl.
Test_x_Lw.interpolate(conditional(le(x1[0],Lw),1.0,0.0))

Indw = [item for item in range(len(Test_x_Lw.dat.data[:])) if Test_x_Lw.dat.data[item] != 0.0]
# yl notes: the set of nodes located before x=Lw. Each element represents the No. of that node.

#print('Update of the solutions:')
PETSc.Sys.Print('Update of the solutions:')
""" *********************************************************************************
    *                                   Time loop                                   *
    ********************************************************************************* """
t_save = t
before_it = time.perf_counter()-start_time # running time from start until this line
t_aux = t
update_wm = 'Yes' # It's used to record the time when t>t_stop for the first time.
smallfac = 10.0**(-10)  #  Onno 13-04-2020

#pdb.set_trace()

while t<=Tend+smallfac: #  while t<Tend-dt: Onno 15-04-2020
    """ *****************************************************************************
        *                               SAVE FUNCTIONS                              *
        ***************************************************************************** """
    if t_save-smallfac < t: # Onno 13-04-2020
        # yl updated:
        progress = format(100*(t-T0)/(Tend-T0), '.3f')+' %'
        tt = format(t, '.3f')
        PETSc.Sys.Print('t= %s, Progress: %s' % (tt, progress))
        #or PETSc.Sys.Print('Progress(%%): ', 100*t/Tend) or PETSc.Sys.Print('Progress', 100*t/Tend,' %%')
        #-------------------------------------------------------------------------------#
        #                                    ENERGY                                     #
        #-------------------------------------------------------------------------------#
        # yl update: XUVW
        energy = save_energy(dim, Lw, H0, g, H, h_n0, psi_1_n0, hat_psi_n0, Uu, Ww, VoW, IoW, WH, b, C11, CN1, CNN, B11, B1N, BN1, BNN, FWF,  
                             A11, AN1, A1N, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN)
        
        if dim=='2D':
            if input_data=='created':
                Energy_file.write('%-25s %-25s %-25s %-25s %-25s %-25s %-25s\n' \
                % (str(t), str(energy), str(h_n0.at(0)), str(psi_1_n0.at(0)), str(hat_psi_n0.at(0)[-1]), str(WM.at(0)), str(dWM_dt.at(0))))
            else:  # input_data=='measurements'
                Energy_file.write('%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s\n' \
                % (str(t), str(energy), str(h_n0.at(10)), str(h_n0.at(20)), str(h_n0.at(40)), str(h_n0.at(49.5)),\
                    str(h_n0.at(50)), str(h_n0.at(54)), str(WM.at(0)), str(dWM_dt.at(0))))
        else: # 3D
            #Energy_file.write('%-25s %-25s %-25s %-25s %-25s %-25s\n' \
            #% (str(t), str(energy), str(h_n0.at(0,0)), str(WM.at(0,0)), str(dWM_dy.at(0,0)), str(dWM_dt.at(0,0))))
            Energy_file.write('%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s\n' \
            % (str(t), str(energy), str(h_n0.at(0,0)), str(h_n0.at(Lx,0)), str(h_n0.at(Lx,Ly)),\
               str(WM.at(0,0)), str(dWM_dy.at(0,0)), str(dWM_dt.at(0,0)), str(WM.at(0,Ly)), str(dWM_dt.at(0,Ly))))
        
        #--output h(x,t) and phi(x,H0,t) and phi(x,0,t) --#
        op_file_name = tt+'.txt'
        op_file = open(os.path.join(save_path, op_file_name), 'w')
        xe = np.linspace(0,Lx,Nx+1)
        for ix in xe:
            op_file.write('%-25s %-25s %-25s %-25s\n' %(str(ix), str(h_n0.at(ix)), str(psi_1_n0.at(ix)), str(hat_psi_n0.at(ix)[-1])))
        op_file.close()

        

        #-------------------------------------------------------------------------------#
        #                               SAVE 3D FUNCTIONS                               #
        #-------------------------------------------------------------------------------#
        if save_pvd:
            #______________________________ Project solutions ______________________________#
            if dim == '2D':
                # To the surface plane (x,y): # add WM by yl.
                x_to_xy(h_n0, WM, psi_1_n0, hat_psi_n0, h_2D, WM_2D, psi_s_2D, psi_i_2D, Indx)
                # In depth (x,y,z):
                for i in range(0,n_z+1):                                     # for each layer
                    phi_projection(i, n_z, waves, Indz, psi_s_2D, psi_i_2D)  # phi(z) = psi_i
                    # WM_3D.dat.data[Indz[i]] = WM.dat.data[0]                 # WM(z) = WM. 
                    # yl: Wrong! will cause mistake in x_transformation
                    # yl updated:
                    WM_3D.dat.data[Indz[i]] = WM_2D.dat.data[:]
            elif dim == '3D':
                # In depth (x,y,z):
                for i in range(0,n_z+1):                                     # for each layer
                    phi_projection(i, n_z, waves, Indz, psi_1_n0, hat_psi_n0)# phi(z) = psi_i
                    WM_3D.dat.data[Indz[i]] = WM.dat.data[:]                     # WM(z) = WM

            #__________________________ Save the fixed coordinates _________________________#
            init_coord = mesh_3D.coordinates.vector().get_local()

            #_________________________________ z-transform _________________________________#
            if dim == '2D':
                z_transform(mesh_3D, n_z, h_2D, beach_s_2D, H0, Indz)
            elif dim == '3D':
                z_transform(mesh_3D, n_z, h_n0, b, H0, Indz)

            #_________________________________ x-transform _________________________________#
            x_transform(mesh_3D, Lw, WM_3D, Indw)

            #_________________________________ Save waves __________________________________#
            save_waves.write(waves)
            # yl notes: now the waves values are saved to the time-dependent mesh.

            #__________________________ Back to the initial mesh ___________________________#
            mesh_3D.coordinates.vector().set_local(init_coord)
        
            #_______________________________ Save wavemaker ________________________________#
            save_WM.write(WM_3D)
            # yl notes: but the wave_maker values are saved to the time-independent mesh.
            
        #_____________________________ Update saving time ______________________________#
        t_save+=dt_save
        

    """ *********************************************************************
        *                            Update time                            *
        ********************************************************************* """

    #_______________________ Update time: t^n -> t^{n+1} _______________________#
    t_half = t+0.5*dt
    t += dt
    
    if input_data=="created":
        if t<=t_stop:                                    # The wavemaker keeps moving
            if scheme=="SV":
                # updated by yl
                WM_expr(WM_half,x,t_half)                            # update R(x,y;t)
                dWM_dt_expr(dWM_half_dt,x,t_half)                    # update dR/dt
                dWM_dy_expr(dWM_half_dy,x,t_half)                    # update dR/dy
                # yl update: XUVW
                Uu_half.assign(Xx*dWM_half_dy)
                Ww_half.assign(Lw-WM_half)
                Vv_half.assign(Lw*Lw+Uu_half*Uu_half)
                VoW_half.assign(Vv_half/Ww_half)
                IoW_half.assign(1/Ww_half)
                WH_half.assign(Ww_half*H)
                XRt_half.assign(Xx*dWM_half_dt)

            # updated by yl
            WM_expr(WM_n1,x,t)                            # update R(x,y;t)
            dWM_dt_expr(dWM_n1_dt,x,t)                    # update dR/dt
            dWM_dy_expr(dWM_n1_dy,x,t)                    # update dR/dy
            # yl update: XUVW
            Ww_n1.assign(Lw-WM_n1)
            
            t_aux = t      # yl notes: store the time when the wave maker stops

        elif t>t_stop and update_wm=='Yes':                       # We stop the wavemaker motion;
            update_wm = 'No'
            if scheme=="SV":
                if t_half<=t_stop:
                    t_aux = t_half
            
            WM_expr(WM_n1,x,t_aux)
            dWM_n1_dt.assign(0.0)
            dWM_dy_expr(dWM_n1_dy,x,t_aux)
            # yl update: XUVW
            Ww_n1.assign(Lw-WM_n1)

            if scheme=="SV":
                WM_expr(WM_half,x,t_aux)
                dWM_half_dt.assign(0.0)
                dWM_dy_expr(dWM_half_dy,x,t_aux)
                # yl update: XUVW
                Uu_half.assign(Xx*dWM_half_dy)
                Ww_half.assign(Lw-WM_half)
                Vv_half.assign(Lw*Lw+Uu_half*Uu_half)
                VoW_half.assign(Vv_half/Ww_half)
                IoW_half.assign(1/Ww_half)
                WH_half.assign(Ww_half*H)
                XRt_half.assign(Xx*dWM_half_dt)
    
    else: # yl added. input_data=='measurements'
        if scheme=="SV":
            WM_expr(WM_half,x,t_half,wm_inter)
            dWM_dt_expr(dWM_half_dt,x,t_half,wm_vel_inter)                    

            # yl update: XUVW
            Ww_half.assign(Lw-WM_half)
            VoW_half.assign(Vv_half/Ww_half)
            IoW_half.assign(1/Ww_half)
            WH_half.assign(Ww_half*H)
            XRt_half.assign(Xx*dWM_half_dt)

        WM_expr(WM_n1,x,t,wm_inter)
        dWM_dt_expr(dWM_n1_dt,x,t,wm_vel_inter)
            
        # yl update: XUVW
        Ww_n1.assign(Lw-WM_n1)

    """ **************************************************
        *            Solve the weak formulations         *
        ************************************************** """
    #___________________ Call the solvers ___________________#
    if scheme=="SE":                     # 1st-order SE scheme
        h_solver.solve()           # get h^{n+1} and hat_psi^*
        psi_solver.solve()                     # get psi^{n+1}
    elif scheme=="SV":                   # 2nd-order SV scheme
        psi_half_solver.solve() # get psi^{n+1/2} and hat_psi^*
        h_solver.solve()        # get h^{n+1} and hat_psi^{**}
        psi_n1_solver.solve()                  # get psi^{n+1}
    
    """ *************************************************
        *               Update the functions            *
        ************************************************* """
    #_________________ Update the solutions ________________#
    h_out, hat_psi_out = w_n1.split()
    h_n0.assign(h_out)
    psi_1_n0.assign(psi_1_n1)
    hat_psi_n0.assign(hat_psi_out) 
    #yl question: why don't update hat_psi?
    # Because we don't need that when we solve for h^n+1 and psi_1^n+1. But we need that when we output .pvd!
    
    #_________________ Update the wavemaker ________________#
    WM.assign(WM_n1)
    dWM_dt.assign(dWM_n1_dt)
    dWM_dy.assign(dWM_n1_dy)

    # yl update: XUVW
    Uu.assign(Xx*dWM_dy)
    Ww.assign(Lw-WM)
    Vv.assign(Lw*Lw+Uu*Uu)
    VoW.assign(Vv/Ww)
    IoW.assign(1/Ww)
    WH.assign(Ww*H)
    XRt.assign(Xx*dWM_dt)

comp_time = time.perf_counter()-start_time
jours = int(comp_time/(24*3600))
heures = int((comp_time-jours*24*3600)/3600)
minutes = int((comp_time-jours*24*3600-heures*3600)/60)
secondes = comp_time -jours*24*3600-heures*3600 - minutes*60
save_README(README_file, Lx, Ly, H0, xb, sb, res_x, Nx, Ny, n_z, gamma, speed, speed, t_stop, Lw, scheme, dt, t,\
            jours, heures, minutes, secondes, comp_time, before_it)
