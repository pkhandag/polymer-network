""" 
Created in: 2021
Purpose: obtain average segment density and total free energy of polymer network with nonlocal inter-segment interactions (using 4-chain model) in 2D
Contact: Pratik Khandagale (pkhandag@andrew.cmu.edu)
"""


#imports
from __future__ import print_function
from fenics import *
from ufl import *
from boxfield import *
from scipy.optimize import fsolve
from numpy.linalg import svd
from sympy import Matrix
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import linalg, matrix
from scipy.integrate import odeint
from tempfile import TemporaryFile
from dolfin import *
from mshr import *
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations_with_replacement  
from numpy import linalg as LA
from scipy.linalg import sqrtm
from xlwt import Workbook 

import numpy as np
import matplotlib.pyplot as plt
import math 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import xlwt 



#####################
## inputs in the model (independent parameters)
#####################

# polymer parameters
N_chain= 100; #number of segments in one single polymer chain in the polymer network
constant_used_in_excluded_volume_parameter= 0.005 # positive for segment repulsion, and vtakes value less than 1
k_B=1.0 # Boltzmann constant (in normalized setting)
Temp=1.0 #temperature of the polymer network  (in normalized setting)

## deformation parameters: principal stretches 
lambda_1= 1
lambda_2= 1
lambda_3= 1

#computational parameters
no_of_elements=35 #no of finite elements alng each X, Y and Z axis 
dt=0.01  #Time step along chain contour. To satisfy CFL numerical stability condition, we need dt < (((x_max_box-x_min_box)/nx_V)**2)/G_chain
c_dirac=0.1 #standard deviation of Gaussian used to approximate Dirac delta function in the initial condition for q and q*


#####################
## other parameters in the model
#####################

a_chain=1/sqrt(N_chain) # segment length
G_chain=((a_chain**2)*N_chain)/6 # defining a constant term ((a_chain**2)*N_chain)/6 in PDE as G_chain
V_seg= a_chain**3 #area of a single segment (area because 2D)
u0=  constant_used_in_excluded_volume_parameter *V_seg ## segment-segment interaction (excluded volume) parameter, This has unit of volume-(to check, look the dirac potential expression)

T=1.0 # final value of chain parameter 's' 
n_t=int(T/dt +1) # number of time steps along the chain contour 
round_const= 12 # no of significant digits  
tol=2e-16 ## tolerance to form submeshes
delta_H_ratio_threshold=1e-3 ## iteration stopping criteria: threshold for relative change in total free energy 
l_RMS= a_chain*(N_chain)**(1/2) ## RMS end-to-end length of one chain


#####################
## initializing variable values
#####################

## initializing total free energies 
W = 0 #total free energy
W_entropic = 0 #entropic free energy
W_interaction = 0 #interaction free energy


## initializing values of relative change in total free energy for checking the stopping criteria for Self Consistent Field Theory iteration 
delta_H_ratio_array=np.zeros(1000) 
delta_H_ratio_array[0]=5
delta_H_ratio=5


############################################################### 
## V mesh forming
############################################################### 

nx_V= no_of_elements #no of finite elements along x-axis
ny_V= no_of_elements #no of finite elements along y-axis
nz_V= no_of_elements #no of finite elements along z-axis

## mesh box size
x_min_box= round( -3* lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )
y_min_box= round( -3* lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )
z_min_box= round( -3* lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )

x_max_box=  round( x_min_box + 6* lambda_1    * ( (1/(sqrt(3)))*l_RMS ), round_const)
y_max_box=  round( y_min_box + 6* lambda_2    * ( (1/(sqrt(3)))*l_RMS ), round_const)
z_max_box=  round( z_min_box + 6* lambda_3    * ( (1/(sqrt(3)))*l_RMS ), round_const)

## domian volume as constant of proportionality for Q and rho 
V_domain= (x_max_box- x_min_box)* (y_max_box- y_min_box)* (z_max_box- z_min_box)


#######################################################


## Define periodic boundary condition
class PeriodicBoundary(SubDomain):

    def inside(self, x, on_boundary):
        # return True if on left or bottom or front boundary or not on the edges
      
        return bool(  (near(x[0], x_min_box) or near(x[1], y_min_box) or near(x[2], z_min_box)) and \
                      ( not    (  (near(x[0], x_min_box) and near(x[1], y_max_box)) or \
                                  (near(x[0], x_min_box) and near(x[2], z_max_box)) or \
                                  (near(x[1], y_min_box) and near(x[0], x_max_box)) or \
                                  (near(x[1], y_min_box) and near(x[2], z_max_box)) or \
                                  (near(x[2], z_min_box) and near(x[0], x_max_box)) or \
                                  (near(x[2], z_min_box) and near(x[1], y_max_box))        )      )    and   on_boundary )
            
            
        # return bool(  (near(x[0], x_min_box) or near(x[1], y_min_box) or near(x[2], z_min_box)) and  on_boundary )
            
            
            
            
    # Map right boundary to left boundary 
    def map(self, x, y):
        if (near(x[0], x_max_box) and near(x[1], y_max_box) and near(x[2], z_max_box)):
            y[0] = x[0] - 2*x_max_box
            y[1] = x[1] - 2*y_max_box
            y[2] = x[2] - 2*z_max_box
        elif (near(x[0], x_max_box) and near(x[1], y_max_box)):
            y[0] = x[0] - 2*x_max_box
            y[1] = x[1] - 2*y_max_box
            y[2] = x[2]
        elif (near(x[1], y_max_box) and near(x[2], z_max_box)):
            y[0] = x[0]
            y[1] = x[1] - 2*y_max_box
            y[2] = x[2] - 2*z_max_box
        elif (near(x[2], z_max_box) and near(x[0], 2*x_max_box)):
            y[0] = x[0] - 2*x_max_box
            y[1] = x[1]
            y[2] = x[2] - 2*z_max_box
        elif near(x[0], x_max_box):
            y[0] = x[0] - 2*x_max_box
            y[1] = x[1]
            y[2] = x[2]
        elif near(x[1], y_max_box):
            y[0] = x[0]
            y[1] = x[1] - 2*y_max_box
            y[2] = x[2]
        elif near(x[2], z_max_box):
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - 2*z_max_box
        else:
            y[0] = 1000.*2*x_max_box
            y[1] = 1000.*2*y_max_box
            y[2] = 1000.*2*z_max_box



## Create mesh and define function space and dof coordinates
mesh = BoxMesh(Point(x_min_box, y_min_box, z_min_box), Point(x_max_box, y_max_box, z_max_box), nx_V, ny_V, nz_V)                                            
V = FunctionSpace(mesh, 'Lagrange', 1, constrained_domain=PeriodicBoundary())        

n_mesh = V.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
d_mesh = mesh.geometry().dim()                                                        
dof_coordinates = V.tabulate_dof_coordinates()                   
dof_coordinates.resize((n_mesh, d_mesh))                                                   
dof_x = dof_coordinates[:, 0]                                                    
dof_y = dof_coordinates[:, 1]     
dof_z = dof_coordinates[:, 2]                                                    

  


################################################################################
## Function for computing q 

def q_computation(X, w):


    
    ##initial q and q_star at t=0
    q_n= Function(V) 
   
    
    ##initial condition for q
    x_cord= X[0]
    y_cord= X[1]
    z_cord= X[2]
    q_0_expression=  Expression( '  ( pow(  (a_chain*sqrt(N_chain)) , 3 ) ) * (1/(sqrt(2*pi)*c_dirac)) * exp(   ( -1/ ( 2*pow(c_dirac,2) ) )  * (    pow(( x[0]- x_cord), 2) +  pow(( x[1]- y_cord), 2) +  pow(( x[2]- z_cord), 2) )   )  '  ,  a_chain= a_chain, N_chain=N_chain, c_dirac=c_dirac, x_cord=x_cord, y_cord=y_cord, z_cord=z_cord, degree=2 )    
    q_0= interpolate(q_0_expression, V)

    ## write initial condition to file
    xdmf_q.write_checkpoint(q_0, "q_label", 0,  XDMFFile.Encoding.HDF5, False)
      
    ##initialize q value at t=0
    q_n.assign(q_0)
    
    
    ######## time stepping for computing q
    for n in range(1,n_t): 
        
        # print(n)
        
        t=dt*n

        
        #defining q, v
        q = TrialFunction(V)
        v = TestFunction(V)
        
        #a and L in fem weak form for fenics (with Crank-Nicolson time stepping)
        a= G_chain*(dt/2)*dot(grad(q), grad(v))*dx + q*v*dx  + dt*w*q*v*dx 
        L= ( q_n*v*dx - G_chain*(dt/2)*dot(grad(q_n), grad(v))*dx  - dt*w*q_n*v*dx)
        
        #solve variational problem
        q = Function(V)
        solve(a == L, q, solver_parameters={'linear_solver': 'gmres', 'preconditioner': 'ilu'})
        
        #saving solution  to file        
        xdmf_q.write_checkpoint(q, "q_label", t,  XDMFFile.Encoding.HDF5, True)
        
        #Update previous solution
        q_n.assign(q)
    
    ######## end of time stepping
    xdmf_q.close()
    
    
    ## returning value of function
    return (q)



        
    
################################################################################
## Function for computing q_star

def q_star_computation(X, w):

    
    ##initial q_star at t=0
    q_star_n= Function(V) 
       
    
    ##initial condition for q_star 
    x_cord= X[0]
    y_cord= X[1]
    z_cord= X[2]
    q_star_0_expression=  Expression( '  ( pow(  (a_chain*sqrt(N_chain)) , 3 ) ) * (1/(sqrt(2*pi)*c_dirac)) * exp(   ( -1/ ( 2*pow(c_dirac,2) ) )  * (    pow(( x[0]- x_cord), 2) +  pow(( x[1]- y_cord), 2) +  pow(( x[2]- z_cord), 2) )   )  '  ,  a_chain= a_chain, N_chain=N_chain, c_dirac=c_dirac, x_cord=x_cord, y_cord=y_cord, z_cord=z_cord, degree=2 )    
    q_star_0= interpolate(q_star_0_expression, V)

    #write
    xdmf_q_star.write_checkpoint(q_star_0, "q_star_label", 0,  XDMFFile.Encoding.HDF5, False)
      
    ##initialize q_star value at t=0
    q_star_n.assign(q_star_0)

    ######## time stepping for q_star
    for n in range(1,n_t): 
        
        # print(n)
        
        t=dt*n
       
        ######################################### computing q_star
        
        #defining q* and v*
        q_star = TrialFunction(V)
        v_star = TestFunction(V)
        
        #a_star and L_star in fem weak form for fenics (with Crank-Nicolson time stepping)
        a_star= G_chain*(dt/2)*dot(grad(q_star), grad(v_star))*dx + q_star*v_star*dx  + dt*w*q_star*v_star*dx 
        L_star= ( q_star_n*v_star*dx - G_chain*(dt/2)*dot(grad(q_star_n), grad(v_star))*dx  - dt*w*q_star_n*v_star*dx)
                
        #solve variational problem
        q_star=Function(V)
        solve(a_star == L_star, q_star, solver_parameters={'linear_solver': 'gmres', 'preconditioner': 'ilu'})
   
        #saving solution to file
        xdmf_q_star.write_checkpoint(q_star, "q_star_label", t,  XDMFFile.Encoding.HDF5, True)
        
         
        #Update previous solution
        q_star_n.assign(q_star)
    
    #### end of time stepping for q_star
    xdmf_q_star.close()

    ## returning values from function
    return (q_star)



################################################################################
   
# Function for single chain computation

def single_chain_computation():
       
    ##computing Q (Complete Partition Function for single chain)
    
    Q=np.zeros(n_t) # Complete Partition Function Q at each position along the chain
    
    phi_chain=Function(V)  # phi function
    phi_chain_temp=Function(V)  # phi function
    phi_chain_numr= phi_chain.vector().get_local() 
    
    for i in range(n_t):
        
        # print(i)

        q_temp = Function(V)
        xdmf_q_call =  XDMFFile("q.xdmf")
        xdmf_q_call.read_checkpoint(q_temp,"q_label",i) 
        xdmf_q_call.close()
        
        q_star_temp = Function(V)
        xdmf_q_star_call =  XDMFFile("q_star.xdmf")
        xdmf_q_star_call.read_checkpoint(q_star_temp,"q_star_label", n_t-1-i) 
        xdmf_q_star_call.close()
       
        Q[i]=assemble((q_temp*q_star_temp)*dx)/V_domain #Q is normalized with dividing by volume of the domain
        
        ## computing average segment density for single chain (phi_chain))  
        q_temp_numr = q_temp.vector().get_local()                                            
        q_star_temp_numr = q_star_temp.vector().get_local()                                            
        phi_chain_temp_numr= phi_chain_temp.vector().get_local() 

        phi_chain_temp_numr= q_temp_numr*q_star_temp_numr
        phi_chain_numr= phi_chain_numr + phi_chain_temp_numr
      
    Q_chain= Q[round(n_t/2)] #Q at s=0.5  
    phi_chain_numr= phi_chain_numr *(1/(V_domain*Q_chain))

    phi_chain.vector().set_local(phi_chain_numr)
    phi_chain.vector().apply('insert')

    ## returning values from function
    return (Q, phi_chain, phi_chain_numr)

################################################################################       



#### computation for initial guess of w


###############################################################################
## 8 Chain computation 
###############################################################################

###############################################################################
## Chain 1 computation
  

########## generating random w

#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))



## Point 1 computation
x1= - round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X1 vector
y1= - round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X1 vector
z1= - round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X1 vector    
X1=np.array([x1,y1,z1]) #position vector of 1st chain start end

# generating vtk files to store q_point_1
# vtkfile_q = File('q_spyder46_3D_solution/q_point_1.pvd')
xdmf_q = XDMFFile("q.xdmf")

q_point_1=Function(V) 
q_point_1= q_computation(X1, w)
 


## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 1st chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0 =Function(V) 
q_point_0 = q_star_computation(X0, w)


Q1, phi_chain_1, phi_chain_1_numr = single_chain_computation()



###############################################################################
## Chain 2 computation
   
   
#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))
   

## Point 2 computation
x2=   round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X2 vector
y2= - round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X2 vector
z2= - round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X2 vector    
X2=np.array([x2,y2,z2]) #position vector of 2nd chain start end

# generating vtk files to store q_point_2
xdmf_q = XDMFFile("q.xdmf")

q_point_2=Function(V) 
q_point_2= q_computation(X2, w)
 

## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 2nd chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0=Function(V) 
q_point_0 = q_star_computation(X0, w)


Q2, phi_chain_2, phi_chain_2_numr = single_chain_computation()


###############################################################################
## Chain 3 computation
   

#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))
   

## Point 3 computation
x3=   round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X3 vector
y3=   round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X3 vector
z3= - round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X3 vector    
X3=np.array([x3,y3,z3]) #position vector of 3rd chain start end

# generating vtk files to store q_point_3
xdmf_q = XDMFFile("q.xdmf")

q_point_3=Function(V) 
q_point_3= q_computation(X3, w)
 

## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 3rd chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0=Function(V) 
q_point_0 = q_star_computation(X0, w)


Q3, phi_chain_3, phi_chain_3_numr = single_chain_computation()


###############################################################################
## Chain 4 computation
   
   
#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))
   

## Point 4 computation
x4= - round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X4 vector
y4=   round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X4 vector
z4= - round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X4 vector    
X4=np.array([x4,y4,z4]) #position vector of 4th chain start end

# generating vtk files to store q_point_4
xdmf_q = XDMFFile("q.xdmf")

q_point_4=Function(V) 
q_point_4= q_computation(X4, w)
 

## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 4th chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0=Function(V) 
q_point_0 = q_star_computation(X0, w)


Q4, phi_chain_4, phi_chain_4_numr = single_chain_computation()


###############################################################################
## Chain 5 computation
   
   
#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))
   

## Point 5 computation
x5= - round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X5 vector
y5= - round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X5 vector
z5=   round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X5 vector    
X5=np.array([x5,y5,z5]) #position vector of 5th chain start end

# generating vtk files to store q_point_1
xdmf_q = XDMFFile("q.xdmf")

q_point_5=Function(V) 
q_point_5= q_computation(X5, w)
 

## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 5th chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0=Function(V) 
q_point_0 = q_star_computation(X0, w)


Q5, phi_chain_5, phi_chain_5_numr = single_chain_computation()


###############################################################################
## Chain 6 computation
   
   
#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))
   

## Point 6 computation
x6=   round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X6 vector
y6= - round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X6 vector
z6=   round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X6 vector    
X6=np.array([x6,y6,z6]) #position vector of 6th chain start end

# generating vtk files to store q_point_6
xdmf_q = XDMFFile("q.xdmf")

q_point_6=Function(V) 
q_point_6= q_computation(X6, w)
 

## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 6th chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0=Function(V) 
q_point_0 = q_star_computation(X0, w)


Q6, phi_chain_6, phi_chain_6_numr = single_chain_computation()


###############################################################################
## Chain 7 computation
   
   
#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))
   

## Point 7 computation
x7=   round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X7 vector
y7=   round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X7 vector
z7=   round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X7 vector    
X7=np.array([x7,y7,z7]) #position vector of 7th chain start end

# generating vtk files to store q_point_7
xdmf_q = XDMFFile("q.xdmf")

q_point_7=Function(V) 
q_point_7= q_computation(X7, w)
 

## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 7th chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0=Function(V) 
q_point_0 = q_star_computation(X0, w)

         
Q7, phi_chain_7, phi_chain_7_numr = single_chain_computation()



###############################################################################
## Chain 8 computation
   
#### generating random w
w = Function(V) 
w.vector().set_local(np.random.random(n_mesh))


## Point 8 computation
x8= - round( lambda_1     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #x component of X8 vector
y8=   round( lambda_2     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #y component of X8 vector
z8=   round( lambda_3     * ( (1/(sqrt(3)))*l_RMS ), round_const )   #z component of X8 vector    
X8=np.array([x8,y8,z8]) #position vector of 8th chain start end

# generating vtk files to store q_point_8
xdmf_q = XDMFFile("q.xdmf")

q_point_8=Function(V) 
q_point_8= q_computation(X8, w)
 

## Point 0 computation
x0= 0 #x component of X0 vector
y0= 0 #y component of X0 vector
z0= 0 #z component of X0 vector
X0=np.array([x0,y0,z0]) #position vector of 8th chain start end

# generating vtk files to store q_point_0
xdmf_q_star = XDMFFile("q_star.xdmf")

q_point_0=Function(V) 
q_point_0 = q_star_computation(X0, w)


Q8, phi_chain_8, phi_chain_8_numr = single_chain_computation()            



###############################################################################
#getting Q at each chain mid point        
  
Q1_mid= Q1[round(n_t/2)] 
Q2_mid= Q2[round(n_t/2)] 
Q3_mid= Q3[round(n_t/2)] 
Q4_mid= Q4[round(n_t/2)]   
Q5_mid= Q5[round(n_t/2)] 
Q6_mid= Q6[round(n_t/2)] 
Q7_mid= Q7[round(n_t/2)] 
Q8_mid= Q8[round(n_t/2)] 



print('Q1 is') 
print(Q1_mid)
print('Q2 is') 
print(Q2_mid)
print('Q3 is') 
print(Q3_mid)
print('Q4 is') 
print(Q4_mid)
print('Q5 is') 
print(Q5_mid)
print('Q6 is') 
print(Q6_mid)
print('Q7 is') 
print(Q7_mid)
print('Q8 is') 
print(Q8_mid)



#computing phi_MF and converting it into fem function

phi_MF=Function(V)  # phi function
phi_MF_numr=phi_MF.vector().get_local()    #phi over mesh nodes(vector of nx.ny.nz size)

phi_MF_numr= phi_chain_1_numr + phi_chain_2_numr + phi_chain_3_numr + phi_chain_4_numr + phi_chain_5_numr +phi_chain_6_numr +phi_chain_7_numr +phi_chain_8_numr 


phi_MF.vector().set_local(phi_MF_numr)
phi_MF.vector().apply('insert') 

  

# free energy
H = k_B*Temp*(1/(2*u0))*assemble(w*w*dx)  - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  +  math.log(Q5_mid)+  math.log(Q6_mid)+  math.log(Q7_mid)+  math.log(Q8_mid)) #total free energy
H_entropic= - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)+  math.log(Q5_mid)+  math.log(Q6_mid)+  math.log(Q7_mid)+  math.log(Q8_mid)  ) #entropic free energy
H_interaction=   k_B*Temp*(1/(2*u0))*assemble(w*w*dx)  #interaction free energy



#############################################################################################
## getting next w
#############################################################################################


# Shift spatially along x the dofs
shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis

class ShiftedExpr_x(UserExpression):
    def __init__(self,func,**kwargs):
        super().__init__(**kwargs)
        self.func = func
    def eval(self,values,x):
        x0_shift = x[0] - shift_length
        if(x0_shift < x_min_box):
            x0_shift += (x_max_box- x_min_box)
        x_shift = np.array([x0_shift, x[1], x[2]])
        values[0] = self.func(x_shift)
    def value_shape(self):
        return ()



class ShiftedExpr_y(UserExpression):
    def __init__(self,func,**kwargs):
        super().__init__(**kwargs)
        self.func = func
    def eval(self,values,x):
        x1_shift = x[1] - shift_length
        if(x1_shift < y_min_box):
            x1_shift += (y_max_box- y_min_box)
        x_shift = np.array([x[0], x1_shift, x[2]])
        values[0] = self.func(x_shift)
    def value_shape(self):
        return ()



class ShiftedExpr_z(UserExpression):
    def __init__(self,func,**kwargs):
        super().__init__(**kwargs)
        self.func = func
    def eval(self,values,x):
        x2_shift = x[2] - shift_length
        if(x2_shift < z_min_box):
            x2_shift += (z_max_box- z_min_box)
        x_shift = np.array([x[0], x[1], x2_shift])
        values[0] = self.func(x_shift)
    def value_shape(self):
        return ()





## middle layer of boxes (in polymer network schematic)
  
phi_MF.set_allow_extrapolation(True)            
shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B13 = interpolate(ShiftedExpr_x(phi_MF),V)
            
shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B15 = interpolate(ShiftedExpr_x(phi_MF),V)
            
phi_MF_B13.set_allow_extrapolation(True)           
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B10 = interpolate(ShiftedExpr_y(phi_MF_B13),V)
          
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B16 = interpolate(ShiftedExpr_y(phi_MF_B13),V)
           
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B11 = interpolate(ShiftedExpr_y(phi_MF),V)
            
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B17 = interpolate(ShiftedExpr_y(phi_MF),V)
           
phi_MF_B15.set_allow_extrapolation(True)           
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B12 = interpolate(ShiftedExpr_y(phi_MF_B15),V)
            
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B18 = interpolate(ShiftedExpr_y(phi_MF_B15),V)
          

## bottom layer of boxes (in polymer network schematic)

shift_length = -2* lambda_3 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B5 = interpolate(ShiftedExpr_z(phi_MF),V)

phi_MF_B5.set_allow_extrapolation(True)            
shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B4 = interpolate(ShiftedExpr_x(phi_MF_B5),V)
         
shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B6 = interpolate(ShiftedExpr_x(phi_MF_B5),V)

phi_MF_B4.set_allow_extrapolation(True)    
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B1 = interpolate(ShiftedExpr_y(phi_MF_B4),V)
   
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B7 = interpolate(ShiftedExpr_y(phi_MF_B4),V)
   
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B2 = interpolate(ShiftedExpr_y(phi_MF_B5),V)
    
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B8 = interpolate(ShiftedExpr_y(phi_MF_B5),V)   

phi_MF_B6.set_allow_extrapolation(True)    
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B3 = interpolate(ShiftedExpr_y(phi_MF_B6),V)
    
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B9 = interpolate(ShiftedExpr_y(phi_MF_B6),V)

    

## top layer of boxes (in polymer network schematic)

shift_length = 2* lambda_3 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B23 = interpolate(ShiftedExpr_z(phi_MF),V)
  
phi_MF_B23.set_allow_extrapolation(True)
shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B22 = interpolate(ShiftedExpr_x(phi_MF_B23),V)
    
shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B24 = interpolate(ShiftedExpr_x(phi_MF_B23),V)   

phi_MF_B22.set_allow_extrapolation(True)  
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B19 = interpolate(ShiftedExpr_y(phi_MF_B22),V)
 
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B25 = interpolate(ShiftedExpr_y(phi_MF_B22),V)
 
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B20 = interpolate(ShiftedExpr_y(phi_MF_B23),V)
      
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B26 = interpolate(ShiftedExpr_y(phi_MF_B23),V)
    
phi_MF_B24.set_allow_extrapolation(True)    
shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B21 = interpolate(ShiftedExpr_y(phi_MF_B24),V)
    
shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
phi_MF_B27 = interpolate(ShiftedExpr_y(phi_MF_B24),V)


## phi_MF_total

phi_MF_total= Function(V)
phi_MF_total_numr=phi_MF_total.vector().get_local()    
phi_MF_total_numr= phi_MF_B1.vector().get_local() + phi_MF_B2.vector().get_local() + phi_MF_B3.vector().get_local() + phi_MF_B4.vector().get_local() + phi_MF_B5.vector().get_local() + phi_MF_B6.vector().get_local() + phi_MF_B7.vector().get_local() + phi_MF_B8.vector().get_local() + phi_MF_B9.vector().get_local() + phi_MF_B10.vector().get_local() + phi_MF_B11.vector().get_local() + phi_MF_B12.vector().get_local() + phi_MF_B13.vector().get_local() + phi_MF.vector().get_local() + phi_MF_B15.vector().get_local() + phi_MF_B16.vector().get_local() + phi_MF_B17.vector().get_local() + phi_MF_B18.vector().get_local() + phi_MF_B19.vector().get_local() + phi_MF_B20.vector().get_local() + phi_MF_B21.vector().get_local() + phi_MF_B22.vector().get_local() + phi_MF_B23.vector().get_local() + phi_MF_B24.vector().get_local() + phi_MF_B25.vector().get_local() + phi_MF_B26.vector().get_local() + phi_MF_B27.vector().get_local()

phi_MF_total.vector().set_local(phi_MF_total_numr)
phi_MF_total.vector().apply('insert') 



# #defining w field for next step
w=Function(V)  
w_numr= w.vector().get_local()

w_numr=  u0* phi_MF_total.vector().get_local()

w.vector().set_local(w_numr)
w.vector().apply('insert')





###############################################################  Iterating for finding equilibrium mean field w
count=0
# while error_w_norm > error_w_norm_threshold:
while abs(delta_H_ratio) > delta_H_ratio_threshold:
    count=count+1
    print(count)
    print(delta_H_ratio)
    

    if count == 51:
        print('count=51')
        # update dt and nt (with the hope of getting convergence less than 1%)
        dt= dt*(0.5)
        n_t=int(T/dt +1)
        
        
    ############################################################################### 
    ## Point 0 computation
  
    # generating vtk files to store q_point_0
    xdmf_q_star = XDMFFile("q_star.xdmf")

    q_point_0=Function(V) 
    q_point_0 = q_star_computation(X0, w)



    ###############################################################################
    ## Chain computation 
    ###############################################################################


    ###############################################################################
    ## Chain 1 computation
    
    # generating vtk files to store q_point_1
    xdmf_q = XDMFFile("q.xdmf")

    q_point_1=Function(V) 
    q_point_1= q_computation(X1, w)
    
    # generating vtk file for phi_chain
    Q1, phi_chain_1, phi_chain_1_numr = single_chain_computation()

    ###############################################################################
    ## Chain 2 computation
    
    # generating vtk files to store q_point_2
    xdmf_q = XDMFFile("q.xdmf")

    q_point_2=Function(V) 
    q_point_2= q_computation(X2, w)
    
    # generating vtk files to store q_point_2
    Q2, phi_chain_2, phi_chain_2_numr= single_chain_computation()
    
    ###############################################################################
    ## Chain 3 computation
    
    # generating vtk files to store q_point_3
    xdmf_q = XDMFFile("q.xdmf")

    q_point_3=Function(V) 
    q_point_3= q_computation(X3, w)
    
    # generating vtk files to store q_point_3
    Q3, phi_chain_3, phi_chain_3_numr= single_chain_computation()        
    
    ###############################################################################
    ## Chain 4 computation
    
    # generating vtk files to store q_point_4
    xdmf_q = XDMFFile("q.xdmf")

    q_point_4=Function(V) 
    q_point_4= q_computation(X4, w)
    
    # generating vtk files to store q_point_4
    Q4, phi_chain_4, phi_chain_4_numr= single_chain_computation()   
    
    
    ###############################################################################
    ## Chain 5 computation
    
    # generating vtk files to store q_point_5
    xdmf_q = XDMFFile("q.xdmf")

    q_point_5=Function(V) 
    q_point_5= q_computation(X5, w)
    
    # generating vtk files to store q_point_5
    Q5, phi_chain_5, phi_chain_5_numr= single_chain_computation()   

    ###############################################################################
    ## Chain 6 computation
    
    # generating vtk files to store q_point_6
    xdmf_q = XDMFFile("q.xdmf")

    q_point_6=Function(V) 
    q_point_6= q_computation(X6, w)
    
    # generating vtk files to store q_point_6
    Q6, phi_chain_6, phi_chain_6_numr= single_chain_computation()   

    ###############################################################################
    ## Chain 7 computation
    
    # generating vtk files to store q_point_7
    xdmf_q = XDMFFile("q.xdmf")

    q_point_7=Function(V) 
    q_point_7= q_computation(X7, w)
    
    # generating vtk files to store q_point_7
    Q7, phi_chain_7, phi_chain_7_numr= single_chain_computation()   

    ###############################################################################
    ## Chain 8 computation
    
    # generating vtk files to store q_point_8
    xdmf_q = XDMFFile("q.xdmf")

    q_point_8=Function(V) 
    q_point_8= q_computation(X8, w)
    
    # generating vtk files to store q_point_8
    Q8, phi_chain_8, phi_chain_8_numr= single_chain_computation()   

            
            
    ###############################################################################
    ##getting Q at each chain mid point        
   
    Q1_mid= Q1[round(n_t/2)] 
    Q2_mid= Q2[round(n_t/2)] 
    Q3_mid= Q3[round(n_t/2)] 
    Q4_mid= Q4[round(n_t/2)] 
    Q5_mid= Q5[round(n_t/2)] 
    Q6_mid= Q6[round(n_t/2)] 
    Q7_mid= Q7[round(n_t/2)] 
    Q8_mid= Q8[round(n_t/2)] 
  
    
    
    print('Q1 is') 
    print(Q1_mid)
    print('Q2 is') 
    print(Q2_mid)
    print('Q3 is') 
    print(Q3_mid)
    print('Q4 is') 
    print(Q4_mid)
    print('Q5 is') 
    print(Q5_mid)
    print('Q6 is') 
    print(Q6_mid)
    print('Q7 is') 
    print(Q7_mid)
    print('Q8 is') 
    print(Q8_mid)




    #computing phi_MF and converting it into fem function

    phi_MF=Function(V)  # phi function
    phi_MF_numr=phi_MF.vector().get_local()    #phi over mesh nodes(vector of nx.ny.nz size)
    phi_MF_numr= phi_chain_1_numr + phi_chain_2_numr + phi_chain_3_numr + phi_chain_4_numr + phi_chain_5_numr +phi_chain_6_numr +phi_chain_7_numr +phi_chain_8_numr 

    phi_MF.vector().set_local(phi_MF_numr)
    phi_MF.vector().apply('insert') 
    
    
    # free energy
    H = k_B*Temp*(1/(2*u0))*assemble(w*w*dx)  - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  +  math.log(Q5_mid)+  math.log(Q6_mid)+  math.log(Q7_mid)+  math.log(Q8_mid)) #total free energy
    H_entropic= - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)+  math.log(Q5_mid)+  math.log(Q6_mid)+  math.log(Q7_mid)+  math.log(Q8_mid)  ) #entropic free energy
    H_interaction=   k_B*Temp*(1/(2*u0))*assemble(w*w*dx) #interaction free energy


    #############################################################################################
    ## getting next w
    #############################################################################################

    
    ## middle layer of boxes
  
    phi_MF.set_allow_extrapolation(True)            
    shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B13 = interpolate(ShiftedExpr_x(phi_MF),V)
                
    shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B15 = interpolate(ShiftedExpr_x(phi_MF),V)
                
    phi_MF_B13.set_allow_extrapolation(True)           
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B10 = interpolate(ShiftedExpr_y(phi_MF_B13),V)
              
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B16 = interpolate(ShiftedExpr_y(phi_MF_B13),V)
               
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B11 = interpolate(ShiftedExpr_y(phi_MF),V)
                
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B17 = interpolate(ShiftedExpr_y(phi_MF),V)
               
    phi_MF_B15.set_allow_extrapolation(True)           
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B12 = interpolate(ShiftedExpr_y(phi_MF_B15),V)
                
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B18 = interpolate(ShiftedExpr_y(phi_MF_B15),V)
              
    
    ## bottom layer of boxes
    
    shift_length = -2* lambda_3 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B5 = interpolate(ShiftedExpr_z(phi_MF),V)
    
    phi_MF_B5.set_allow_extrapolation(True)            
    shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B4 = interpolate(ShiftedExpr_x(phi_MF_B5),V)
             
    shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B6 = interpolate(ShiftedExpr_x(phi_MF_B5),V)
    
    phi_MF_B4.set_allow_extrapolation(True)    
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B1 = interpolate(ShiftedExpr_y(phi_MF_B4),V)
       
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B7 = interpolate(ShiftedExpr_y(phi_MF_B4),V)
       
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B2 = interpolate(ShiftedExpr_y(phi_MF_B5),V)
        
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B8 = interpolate(ShiftedExpr_y(phi_MF_B5),V)   
    
    phi_MF_B6.set_allow_extrapolation(True)    
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B3 = interpolate(ShiftedExpr_y(phi_MF_B6),V)
        
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B9 = interpolate(ShiftedExpr_y(phi_MF_B6),V)
    
        
    
    ## top layer of boxes
    
    shift_length = 2* lambda_3 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B23 = interpolate(ShiftedExpr_z(phi_MF),V)
      
    phi_MF_B23.set_allow_extrapolation(True)
    shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B22 = interpolate(ShiftedExpr_x(phi_MF_B23),V)
        
    shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B24 = interpolate(ShiftedExpr_x(phi_MF_B23),V)   
    
    phi_MF_B22.set_allow_extrapolation(True)  
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B19 = interpolate(ShiftedExpr_y(phi_MF_B22),V)
     
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B25 = interpolate(ShiftedExpr_y(phi_MF_B22),V)
 
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B20 = interpolate(ShiftedExpr_y(phi_MF_B23),V)
          
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B26 = interpolate(ShiftedExpr_y(phi_MF_B23),V)
        
    phi_MF_B24.set_allow_extrapolation(True)    
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B21 = interpolate(ShiftedExpr_y(phi_MF_B24),V)
        
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B27 = interpolate(ShiftedExpr_y(phi_MF_B24),V)
    

    ## phi_MF_total
    
    phi_MF_total= Function(V)
    phi_MF_total_numr=phi_MF_total.vector().get_local()    
    phi_MF_total_numr= phi_MF_B1.vector().get_local() + phi_MF_B2.vector().get_local() + phi_MF_B3.vector().get_local() + phi_MF_B4.vector().get_local() + phi_MF_B5.vector().get_local() + phi_MF_B6.vector().get_local() + phi_MF_B7.vector().get_local() + phi_MF_B8.vector().get_local() + phi_MF_B9.vector().get_local() + phi_MF_B10.vector().get_local() + phi_MF_B11.vector().get_local() + phi_MF_B12.vector().get_local() + phi_MF_B13.vector().get_local() + phi_MF.vector().get_local() + phi_MF_B15.vector().get_local() + phi_MF_B16.vector().get_local() + phi_MF_B17.vector().get_local() + phi_MF_B18.vector().get_local() + phi_MF_B19.vector().get_local() + phi_MF_B20.vector().get_local() + phi_MF_B21.vector().get_local() + phi_MF_B22.vector().get_local() + phi_MF_B23.vector().get_local() + phi_MF_B24.vector().get_local() + phi_MF_B25.vector().get_local() + phi_MF_B26.vector().get_local() + phi_MF_B27.vector().get_local()
    
    phi_MF_total.vector().set_local(phi_MF_total_numr)
    phi_MF_total.vector().apply('insert') 

    
    ############################################################### 
    ## delta_H_ratio check
    ############################################################### 
    
    if count != 1 and abs((H-H_next)/H_next) < delta_H_ratio_threshold:
        H_next=H
        H_next_entropic= H_entropic
        H_next_interaction= H_interaction
        print('iteration ended after break')
        break
    
    
    
    ############################################################### 
    ## getting w for next time step
    ###############################################################                 
    
    #defining w field for next step
    w_next=Function(V)  # phi function
    w_next_numr= w_next.vector().get_local()
    
    w_next_numr=  u0* phi_MF_total.vector().get_local()
    
    w_next.vector().set_local(w_next_numr)
    w_next.vector().apply('insert')
      
    
    ###############################################################################
    ## computation for w_next
    ###############################################################################


    ############################################################################### 
    ## Point 0 computation
  
    # generating vtk files to store q_point_0
    xdmf_q_star = XDMFFile("q_star.xdmf")

    q_point_0=Function(V) 
    q_point_0 = q_star_computation(X0, w_next)



    ###############################################################################
    ## Chain computation 
    ###############################################################################


    ###############################################################################
    ## Chain 1 computation
    
    # generating vtk files to store q_point_1
    xdmf_q = XDMFFile("q.xdmf")

    q_point_1=Function(V) 
    q_point_1= q_computation(X1, w_next)
    
    # generating vtk file for phi_chain
    Q1, phi_chain_1, phi_chain_1_numr = single_chain_computation()

    ###############################################################################
    ## Chain 2 computation
    
    # generating vtk files to store q_point_2
    xdmf_q = XDMFFile("q.xdmf")

    q_point_2=Function(V) 
    q_point_2= q_computation(X2, w_next)
    
    # generating vtk files to store q_point_2
    Q2, phi_chain_2, phi_chain_2_numr= single_chain_computation()
    
    ###############################################################################
    ## Chain 3 computation
    
    # generating vtk files to store q_point_3
    xdmf_q = XDMFFile("q.xdmf")

    q_point_3=Function(V) 
    q_point_3= q_computation(X3, w_next)
    
    # generating vtk files to store q_point_3
    Q3, phi_chain_3, phi_chain_3_numr= single_chain_computation()        
    
    ###############################################################################
    ## Chain 4 computation
    
    # generating vtk files to store q_point_4
    xdmf_q = XDMFFile("q.xdmf")

    q_point_4=Function(V) 
    q_point_4= q_computation(X4, w_next)
    
    # generating vtk files to store q_point_4
    Q4, phi_chain_4, phi_chain_4_numr= single_chain_computation()   
    
    
    ###############################################################################
    ## Chain 5 computation
    
    # generating vtk files to store q_point_5
    xdmf_q = XDMFFile("q.xdmf")

    q_point_5=Function(V) 
    q_point_5= q_computation(X5, w_next)
    
    # generating vtk files to store q_point_5
    Q5, phi_chain_5, phi_chain_5_numr= single_chain_computation()   

    ###############################################################################
    ## Chain 6 computation
    
    # generating vtk files to store q_point_6
    xdmf_q = XDMFFile("q.xdmf")

    q_point_6=Function(V) 
    q_point_6= q_computation(X6, w_next)
    
    # generating vtk files to store q_point_6
    Q6, phi_chain_6, phi_chain_6_numr= single_chain_computation()   

    ###############################################################################
    ## Chain 7 computation
    
    # generating vtk files to store q_point_7
    xdmf_q = XDMFFile("q.xdmf")

    q_point_7=Function(V) 
    q_point_7= q_computation(X7, w_next)
    
    # generating vtk files to store q_point_7
    Q7, phi_chain_7, phi_chain_7_numr= single_chain_computation()   

    ###############################################################################
    ## Chain 8 computation
    
    # generating vtk files to store q_point_8
    xdmf_q = XDMFFile("q.xdmf")

    q_point_8=Function(V) 
    q_point_8= q_computation(X8, w_next)
    
    # generating vtk files to store q_point_8
    Q8, phi_chain_8, phi_chain_8_numr= single_chain_computation()   

 
    ###############################################################################
    ##getting Q at each chain mid point        
   
    Q1_mid= Q1[round(n_t/2)] 
    Q2_mid= Q2[round(n_t/2)] 
    Q3_mid= Q3[round(n_t/2)] 
    Q4_mid= Q4[round(n_t/2)] 
    Q5_mid= Q5[round(n_t/2)] 
    Q6_mid= Q6[round(n_t/2)] 
    Q7_mid= Q7[round(n_t/2)] 
    Q8_mid= Q8[round(n_t/2)] 
  
  
    
    print('Q1 is') 
    print(Q1_mid)
    print('Q2 is') 
    print(Q2_mid)
    print('Q3 is') 
    print(Q3_mid)
    print('Q4 is') 
    print(Q4_mid)
    print('Q5 is') 
    print(Q5_mid)
    print('Q6 is') 
    print(Q6_mid)
    print('Q7 is') 
    print(Q7_mid)
    print('Q8 is') 
    print(Q8_mid)



    #computing phi_MF and converting it into fem function

    phi_MF=Function(V)  # phi function
    phi_MF_numr=phi_MF.vector().get_local()    #phi over mesh nodes(vector of nx.ny.nz size)
    phi_MF_numr= phi_chain_1_numr + phi_chain_2_numr + phi_chain_3_numr + phi_chain_4_numr + phi_chain_5_numr +phi_chain_6_numr +phi_chain_7_numr +phi_chain_8_numr 


    phi_MF.vector().set_local(phi_MF_numr)
    phi_MF.vector().apply('insert') 
    
    # free energy
    H_next = k_B*Temp*(1/(2*u0))*assemble(w_next*w_next*dx)  - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  +  math.log(Q5_mid)+  math.log(Q6_mid)+  math.log(Q7_mid)+  math.log(Q8_mid))   # total free energy    
    H_next_entropic= - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)+  math.log(Q5_mid)+  math.log(Q6_mid)+  math.log(Q7_mid)+  math.log(Q8_mid)  ) #entropic free energy
    H_next_interaction=   k_B*Temp*(1/(2*u0))*assemble(w_next*w_next*dx) #interaction free energy



    # computing relative change in H and H_next
    delta_H_ratio=(H_next-H)/H
    delta_H_ratio_array[count]=delta_H_ratio
    
    print(delta_H_ratio)

    if abs((H_next-H)/H) < delta_H_ratio_threshold:
        print('iteration ended after break')
        break
    


    #############################################################################################
    ## getting next w
    #############################################################################################

    
    ## middle layer of boxes
  
    phi_MF.set_allow_extrapolation(True)            
    shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B13 = interpolate(ShiftedExpr_x(phi_MF),V)
                
    shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B15 = interpolate(ShiftedExpr_x(phi_MF),V)
                
    phi_MF_B13.set_allow_extrapolation(True)           
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B10 = interpolate(ShiftedExpr_y(phi_MF_B13),V)
              
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B16 = interpolate(ShiftedExpr_y(phi_MF_B13),V)
               
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B11 = interpolate(ShiftedExpr_y(phi_MF),V)
                
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B17 = interpolate(ShiftedExpr_y(phi_MF),V)
               
    phi_MF_B15.set_allow_extrapolation(True)           
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B12 = interpolate(ShiftedExpr_y(phi_MF_B15),V)
                
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B18 = interpolate(ShiftedExpr_y(phi_MF_B15),V)
              
    
    ## bottom layer of boxes
    
    shift_length = -2* lambda_3 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B5 = interpolate(ShiftedExpr_z(phi_MF),V)
    
    phi_MF_B5.set_allow_extrapolation(True)            
    shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B4 = interpolate(ShiftedExpr_x(phi_MF_B5),V)
             
    shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B6 = interpolate(ShiftedExpr_x(phi_MF_B5),V)
    
    phi_MF_B4.set_allow_extrapolation(True)    
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B1 = interpolate(ShiftedExpr_y(phi_MF_B4),V)
       
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B7 = interpolate(ShiftedExpr_y(phi_MF_B4),V)
       
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B2 = interpolate(ShiftedExpr_y(phi_MF_B5),V)
        
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B8 = interpolate(ShiftedExpr_y(phi_MF_B5),V)   
    
    phi_MF_B6.set_allow_extrapolation(True)    
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B3 = interpolate(ShiftedExpr_y(phi_MF_B6),V)
        
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B9 = interpolate(ShiftedExpr_y(phi_MF_B6),V)
    
        
    
    ## top layer of boxes
    
    shift_length = 2* lambda_3 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B23 = interpolate(ShiftedExpr_z(phi_MF),V)
      
    phi_MF_B23.set_allow_extrapolation(True)
    shift_length = -2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B22 = interpolate(ShiftedExpr_x(phi_MF_B23),V)
        
    shift_length = 2* lambda_1 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B24 = interpolate(ShiftedExpr_x(phi_MF_B23),V)   
    
    phi_MF_B22.set_allow_extrapolation(True)  
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B19 = interpolate(ShiftedExpr_y(phi_MF_B22),V)
     
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B25 = interpolate(ShiftedExpr_y(phi_MF_B22),V)
 
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B20 = interpolate(ShiftedExpr_y(phi_MF_B23),V)
          
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B26 = interpolate(ShiftedExpr_y(phi_MF_B23),V)
        
    phi_MF_B24.set_allow_extrapolation(True)    
    shift_length = -2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B21 = interpolate(ShiftedExpr_y(phi_MF_B24),V)
        
    shift_length = 2* lambda_2 * ( (1/(sqrt(3)))*l_RMS )  #shift along x-axis
    phi_MF_B27 = interpolate(ShiftedExpr_y(phi_MF_B24),V)
    

    ## phi_MF_total
    
    phi_MF_total= Function(V)
    phi_MF_total_numr=phi_MF_total.vector().get_local()    
    phi_MF_total_numr= phi_MF_B1.vector().get_local() + phi_MF_B2.vector().get_local() + phi_MF_B3.vector().get_local() + phi_MF_B4.vector().get_local() + phi_MF_B5.vector().get_local() + phi_MF_B6.vector().get_local() + phi_MF_B7.vector().get_local() + phi_MF_B8.vector().get_local() + phi_MF_B9.vector().get_local() + phi_MF_B10.vector().get_local() + phi_MF_B11.vector().get_local() + phi_MF_B12.vector().get_local() + phi_MF_B13.vector().get_local() + phi_MF.vector().get_local() + phi_MF_B15.vector().get_local() + phi_MF_B16.vector().get_local() + phi_MF_B17.vector().get_local() + phi_MF_B18.vector().get_local() + phi_MF_B19.vector().get_local() + phi_MF_B20.vector().get_local() + phi_MF_B21.vector().get_local() + phi_MF_B22.vector().get_local() + phi_MF_B23.vector().get_local() + phi_MF_B24.vector().get_local() + phi_MF_B25.vector().get_local() + phi_MF_B26.vector().get_local() + phi_MF_B27.vector().get_local()
    
    phi_MF_total.vector().set_local(phi_MF_total_numr)
    phi_MF_total.vector().apply('insert') 

 
    
    ############################################################### 
    ## getting w for next time step
    ###############################################################                 
    
    #defining w field for next step
    w=Function(V)  # phi function
    w_numr= w.vector().get_local()
    
    w_numr=  u0* phi_MF_total.vector().get_local()
    
    w.vector().set_local(w_numr)
    w.vector().apply('insert')

  



## converged free energy 
W= H_next # total free energy
W_entropic= H_next_entropic # entropic free energy
W_interaction= H_next_interaction #interaction free energy



############################ End of computation




##############################################3
## saving avg segment density results
##############################################3


File('phi_MF_B1.pvd') << (phi_MF_B1)
File('phi_MF_B2.pvd') << (phi_MF_B2)
File('phi_MF_B3.pvd') << (phi_MF_B3)
File('phi_MF_B4.pvd') << (phi_MF_B4)
File('phi_MF_B5.pvd') << (phi_MF_B5)
File('phi_MF_B6.pvd') << (phi_MF_B6)
File('phi_MF_B7.pvd') << (phi_MF_B7)
File('phi_MF_B8.pvd') << (phi_MF_B8)
File('phi_MF_B9.pvd') << (phi_MF_B9)
File('phi_MF_B10.pvd') << (phi_MF_B10)
File('phi_MF_B11.pvd') << (phi_MF_B11)
File('phi_MF_B12.pvd') << (phi_MF_B12)
File('phi_MF_B13.pvd') << (phi_MF_B13)
File('phi_MF.pvd') << (phi_MF)
File('phi_MF_B15.pvd') << (phi_MF_B15)
File('phi_MF_B16.pvd') << (phi_MF_B16)
File('phi_MF_B17.pvd') << (phi_MF_B17)
File('phi_MF_B18.pvd') << (phi_MF_B18)
File('phi_MF_B19.pvd') << (phi_MF_B19)
File('phi_MF_B20.pvd') << (phi_MF_B20)
File('phi_MF_B21.pvd') << (phi_MF_B21)
File('phi_MF_B22.pvd') << (phi_MF_B22)
File('phi_MF_B23.pvd') << (phi_MF_B23)
File('phi_MF_B24.pvd') << (phi_MF_B24)
File('phi_MF_B25.pvd') << (phi_MF_B25)
File('phi_MF_B26.pvd') << (phi_MF_B26)
File('phi_MF_B27.pvd') << (phi_MF_B27)
File('phi_MF_total.pvd') << (phi_MF_total)

File('w.pvd') << (w)
File('w_next.pvd') << (w_next)






