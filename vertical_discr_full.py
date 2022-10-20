# full version

from sympy import *

""" *************************************************************
    *                      Lagrange polynomial                  *
    *************************************************************
    This function gives the expression of the Lagrange polynomial
    varphi_i(z) of order n, for z between z=0 and z=H_0.         """

def varphi_expr(i, n, H0):
    z=Symbol('z')                                     # z-coordinate
    k = Symbol('k')                         # index k in the product
    z_k = H0*((n-k)/n)                    # discrete coordinates z_k
    sigma = lambdify(k,z_k,"numpy")                 # sigma(k) = z_k
    varphi_z = \
    (Product((z-sigma(k))/(sigma(i)-sigma(k)),(k,0,i-1))\
     *Product((z-sigma(k))/(sigma(i)-sigma(k)),(k,i+1,n))).doit()
    return varphi_z

""" ************************************
    *           d(\varphi)/dz          *
    ************************************
    This function returns the expression
    of the z--derivative of the Lagrange
    polynomial, that is d(\varphi)/dz . """

def deriv_varphi_expr(varphi_expr):
    z = Symbol('z')       # coordinate z
    deriv = diff(varphi_expr,z)
    return deriv         # d(\varphi(z))/dz

""" ***********************************************
    *               Vertical matrices             *
    *********************************************** """

def M_ij(i,j,n,H0):
    z=Symbol('z')
    expr_M = varphi_expr(i,n,H0)*varphi_expr(j,n,H0)
    M = integrate(expr_M, (z,0,H0))
    return M

def A_ij(i,j,n,H0):
    z=Symbol('z')
    expr_A = deriv_varphi_expr(varphi_expr(i,n,H0))\
        *deriv_varphi_expr(varphi_expr(j,n,H0))
    A = integrate(expr_A, (z,0,H0))
    return A

# yl added B_ij: for full weak forms
def B_ij(i,j,n,H0):
    z=Symbol('z')
    expr_B = varphi_expr(i,n,H0)*deriv_varphi_expr(varphi_expr(j,n,H0))
    B = integrate(expr_B, (z,0,H0))
    return B

def D_ij(i,j,n,H0):
    z=Symbol('z')
    expr_D = z*varphi_expr(i,n,H0)\
        *deriv_varphi_expr(varphi_expr(j,n,H0))
    D = integrate(expr_D, (z,0,H0))
    return D

# yl added C_ij: for full weak forms
def C_ij(i,j,n,H0):
    z=Symbol('z')
    expr_C = z*deriv_varphi_expr(varphi_expr(i,n,H0))\
        *deriv_varphi_expr(varphi_expr(j,n,H0))
    C = integrate(expr_C, (z,0,H0))
    return C

def S_ij(i,j,n,H0):
    z=Symbol('z')
    expr_S = z*z*deriv_varphi_expr(varphi_expr(i,n,H0))\
        *deriv_varphi_expr(varphi_expr(j,n,H0))
    S = integrate(expr_S, (z,0,H0))
    return S

def I_i(i,n,H0):
    z=Symbol('z')
    expr_I = varphi_expr(i,n,H0)
    I = integrate(expr_I, (z,0,H0))
    return I