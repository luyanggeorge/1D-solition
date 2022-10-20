from firedrake import *
import numpy as np

"""
    *********************************************
    *                 Test case                 *
    *********************************************"""
def test_case():
    #________________ Kind of data ________________#
    #input_data = "measurements"  # from experiments
    input_data = "created"       # set the wavemaker
    #______________ Temporal scheme _______________#
    scheme = "SE"
    #  "SE": Symplectic-Euler ; "SV": Stormer-Verlet
    #__________________ Dimension _________________#
    dim = "2D"
    #"2D": R(t) and b(x); "3D": R(y,t) and/or b(x,y)
    # if input = measurements, the dim must be 2D.
    #______ Path and name of the saved files ______#
    save_path = 'data/'+scheme+'/'+dim+'/KPE_sech_2Aug/'
    # ----yl added. whether the seabed is flat or not
    bottom = 'flat' 
    # 'flat':b(x,y)=0; 'nonuniform':b(x,y)!=0
    # ----yl added. Whether or not to apply mild-slope approximation (MSA)
    FWF = 0
    # 1: use full weak forms (FWF); 0: use mild-slope approximations.
    save_pvd = False
    # Whether or not to save the 3D results into pvd files 
    return input_data, scheme, dim, save_path, bottom, FWF, save_pvd

# loading data from files. yl added.
def load_wavemaker(dt):
    wm_data_0     = np.loadtxt('202002/PistonMotion.dat',usecols=1) # measured motion
    t_data_0      = np.loadtxt('202002/PistonVelocity.dat',usecols=0) # measured time
    wm_vel_data_0 = np.loadtxt('202002/PistonVelocity.dat',usecols=1) # measured velocity
    # add an value to avoid overflow in the time loop
    t_data      = np.append(t_data_0,[t_data_0[-1]+2*dt]) 
    wm_data     = np.append(wm_data_0,[wm_data_0[-1]])
    wm_vel_data = np.append(wm_vel_data_0,[wm_vel_data_0[-1]])
    return t_data, wm_data, wm_vel_data

"""
    *********************************************************************
    *                         Numerical domain                          *
    *********************************************************************"""
def domain(bottom):
    #______________________ Beach ______________________#
    H0 = 20.0
    xb = 0.0                                          # Start of the beach
    sb = 0.0                                           # Slope of the beach
    # yl update
    def H_expr(V,x):
        return Function(V).interpolate(H0-conditional(le(x[0],xb),0.0,sb*(x[0]-xb)))
    #______________________ Basin ______________________#
    if bottom=='nonuniform':
        Hend = 0.5                              # Depth at the end of the beach
        Lx = xb +(H0-Hend)/sb                                 # Length in x
    else:
        Hend = H0
        Lx = 8000
    Ly = 1.0                                                  # Length in y
    Lw = 1.0                                        # End of the x-transform
    res_x = 0.5                                             # x-resolution
    res_y = 1.0                                            # y-resolution
    n_z = 5                                         # Order of the expansion
    return H0, xb, sb, H_expr, Hend, Lx, Ly, Lw, res_x, res_y, n_z


"""
    **************************************************************************
    *                                Wavemaker                               *
    **************************************************************************"""
def wavemaker(dim, H0, Ly, Lw, input_data):
    #_____________________________ Characteristics _____________________________#
    g = 9.81                                             # Gravitational constant
    gamma = 0.0                                                  # Wave amplitude                                     
    w = 1.0

    # yl update
    if input_data=='created':
        if dim == "2D":
            def WM_expr(function,x,t):
                function.interpolate(conditional(le(x[0],Lw),-gamma*cos(w*t),0.0))
            def dWM_dt_expr(function,x,t):
                function.interpolate(conditional(le(x[0],Lw),gamma*w*sin(w*t),0.0))
            def dWM_dy_expr(function,x,t):
                function.assign(0.0)
        elif dim == "3D":
            def WM_expr(function,x,t):
                function.interpolate(conditional(le(x[0],Lw),-gamma*((x[1]-0.5*Ly)/(0.5*Ly))*cos(w*t),0.0))
            def dWM_dt_expr(function,x,t):
                function.interpolate(conditional(le(x[0],Lw),gamma*w*((x[1]-0.5*Ly)/(0.5*Ly))*sin(w*t),0.0))
            def dWM_dy_expr(function,x,t):
                function.interpolate(conditional(le(x[0],Lw),-gamma*cos(w*t)/(0.5*Ly),0.0))
    
    else: # interpolations of experimental data
        def WM_expr(function, x, t, inter_data): 
            R_t = round(float(inter_data(t)),4)
            function.interpolate(conditional(le(x[0],Lw),R_t,0.0))
        def dWM_dt_expr(function, x, t, inter_data):
            Rt_t= round(float(inter_data(t)),4)
            function.interpolate(conditional(le(x[0],Lw),Rt_t,0.0))
        def dWM_dy_expr(function):
            function.assign(0.0)

    # initialise the system with 2D exact solution of KPE
    A  = 0.1
    e  = 0.05
    mu = e*e
    a_eta = A*H0*e
    a_phi = e*H0*np.sqrt(g*H0/mu)
    speed = (1+0.5*A*e)*sqrt(g*H0)

    def xi_expr(function,x,t):
        function.interpolate( (sqrt(3*e*A)/(2*H0))*(x[0] - speed*t) )

    def h_ex_expr(function,xi): # h = H0 + eta
        function.assign( H0 + a_eta*(pow(cosh(xi),-2)) )

    def Phi_expr(function,xi):
        function.assign( sqrt(4*e*A/3)*(tanh(xi)+1) )

    def Phi_2x_expr(function,xi):
        function.assign( -(A*e/mu) * sqrt(3*e*A) * tanh(xi) * pow(cosh(xi),-2) )

    def Phi_4x_expr(function,xi):
        function.assign( 3 * pow(A*e/mu,2) * sqrt(3*e*A) *\
                             (2*tanh(xi)*pow(cosh(xi),-4)-pow(tanh(xi),3)*pow(cosh(xi),-2)) )

    def phi_ex_expr(function,Phi,Phi_2x,Phi_4x,z):
        function.assign( a_phi*(Phi - 0.5*mu*pow(z/H0,2)*Phi_2x + (mu*mu/24)*pow(z/H0,4)*Phi_4x) )

    return g, gamma, speed, WM_expr, dWM_dt_expr, dWM_dy_expr,\
           xi_expr, h_ex_expr, Phi_expr, Phi_2x_expr, Phi_4x_expr, phi_ex_expr


"""
    ***********************************
    *               Time              *
    ***********************************"""
def set_time(speed, Lx, res_x, Ly, res_y, H0, n_z):
    Nx = round(Lx/res_x)    # Number of elements in x (round to the nearest integer)
    dx = Lx/Nx
    dt = 0.02                      # time step
    T0 = round(0.5*Lx/speed)       # Initial time
    Tend = T0 + 10 #0.25*Lx/speed      # Final time 
    t_stop = Tend
    t = T0                  # Temporal variable 
    dt_save = 5*dt          # saving time step
    return T0, t, dt, Tend, t_stop, dt_save