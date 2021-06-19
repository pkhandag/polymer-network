""" 
Created in: 2021
Purpose: evaluate average segment density and total free energy of 2D polymer network with nonlocal inter-segment interactions 
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
## inputs in the model(independent parameters)
#####################

#polymer parameters
N_chain= 100; #number of segments in one single polymer chain in the polymer network
constant_used_in_excluded_volume_parameter= 0.1 # positive for segment repulsion, and vtakes value less than 1
k_B=1.0 # Boltzmann constant (in normalized setting)
Temp=1.0 #temperature of the polymer network  (in normalized setting)

## no of deformed states for polymer network considered in the code
no_of_stretch_values=1 

# computation parameters
dt=0.01  #Time step along chain contour. CFL numerical stability condition: dt < (((x_max_box-x_min_box)/nx_V)**2)/G_chain
no_elements_per_lengthscale=20 #no of finite elements per lengthscale (1 lengthscale= a*N^(1/2))
c_dirac=0.1 #standard deviation of Gaussian used to approximate Dirac delta function in the initial condition for q and q*

#####################
## other parameters in the model
#####################

a_chain=1/sqrt(N_chain) # segment length
G_chain=((a_chain**2)*N_chain)/6 # defining constant ((a_chain**2)*N_chain)/6 in PDE as G_chain
V_seg= a_chain**2 #area of a single segment (area because 2D)
u0=  constant_used_in_excluded_volume_parameter *V_seg ## segment-segment interaction (excluded volume) parameter, This has unit of volume-(to check, look the dirac potential expression)

T=1.0 # final value of chain parameter 's' 
n_t=int(T/dt +1) # number of time steps along the chain contour 
round_const= 12 # no of significant digits  
tol=2e-16 ## tolerance to form submeshes
delta_H_ratio_threshold=1e-3 ## iteration stopping criteria: threshold for relative change in total free energy 

#####################
## initializing variable values
#####################

## initializing free energy values 
W = np.zeros(shape=(no_of_stretch_values, no_of_stretch_values)) #total free energy (entropic+interaction)
W_entropic = np.zeros(shape=(no_of_stretch_values, no_of_stretch_values)) #entropic free energy
W_interaction = np.zeros(shape=(no_of_stretch_values, no_of_stretch_values)) #interaction free energy

## initializing principal stretches
lambda_1_stretch=np.zeros(no_of_stretch_values)
lambda_2_stretch=np.zeros(no_of_stretch_values)

## initializing values of relative change in free energy for checking Self Consistent Field Theory iteration stopping criteria
delta_H_ratio_array=np.zeros(1000) #initial delta_H_ratio 
delta_H_ratio_array[0]=5
delta_H_ratio=5

##################################### 


for i1 in range(no_of_stretch_values):
    for i2 in range(no_of_stretch_values):
           
        
        ## principal stretch values
        lambda_1_stretch[i1]= 1
        lambda_2_stretch[i2]= 1      
    
 
        ############################################################### 
        ## V mesh forming
        ############################################################### 
    
        # no of finite elements
        nx_V =   3* int( round( no_elements_per_lengthscale* (lambda_1_stretch[i1])      *2*  (  1/(2*sqrt(2))   ) ) )    #no of elements along x # 42 for lambda_1=1
        ny_V =   3* int( round( no_elements_per_lengthscale* (lambda_2_stretch[i2])      *2*  (  1/(2*sqrt(2))   ) ) )    #no of elements along y
        
        # nx_V= 90
        # ny_V= 90
        
        ## mesh box size
        x_min_box= round( -3* lambda_1_stretch[i1]     * ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_min_box= round( -3* lambda_2_stretch[i2] * ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
    
        x_max_box=  round( x_min_box + 6* lambda_1_stretch[i1]     * ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box=  round( y_min_box + 6* lambda_2_stretch[i2]     * ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        ## domian volume as constant of proportionality for Q and rho 
        V_domain= (x_max_box- x_min_box)* (y_max_box- y_min_box)
        

                
      

        #######################################################

        

        
        ## Define periodic boundary condition
        class PeriodicBoundary(SubDomain):
        
            def inside(self, x, on_boundary):
                # return True if on left or bottom boundary
                # return bool((near(x[0], x_min_box) or near(x[1], y_min_box) ) and on_boundary)
    
                return bool( (near(x[0], x_min_box) or near(x[1], y_min_box) )          and \
                             ( not(    (near(x[0], x_min_box) and near(x[1], y_max_box)) or \
                                       (near(x[0], x_max_box) and near(x[1], y_min_box))       )     )  and on_boundary  )
            
            def map(self, x, y):
    
                if near(x[0], x_max_box) and near(x[1], y_max_box):
                    y[0] = x[0] - 2*x_max_box
                    y[1] = x[1] - 2*y_max_box
                elif near(x[0], x_max_box):
                    y[0] = x[0] - 2*x_max_box
                    y[1] = x[1]
                elif near(x[1], y_max_box):
                    y[0] = x[0]
                    y[1] = x[1] - 2*y_max_box
                else:
                    y[0] = 1000*2*x_max_box
                    y[1] = 1000*2*y_max_box 
        
        
        ## Create mesh and define function space and dof coordinates
        mesh = RectangleMesh(Point(x_min_box, y_min_box), Point(x_max_box, y_max_box), nx_V, ny_V)                                            
        V = FunctionSpace(mesh, 'Lagrange', 1, constrained_domain=PeriodicBoundary())        
        
        n_mesh = V.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh = mesh.geometry().dim()                                                        
        dof_coordinates = V.tabulate_dof_coordinates()                   
        dof_coordinates.resize((n_mesh, d_mesh))                                                   
        dof_x = dof_coordinates[:, 0]                                                    
        dof_y = dof_coordinates[:, 1]     
        
        
        
        
        ############################################################### 
        ## B13 mesh forming
        ############################################################### 
            
        nx_B13= nx_V
        ny_B13= ny_V
        
        ## mesh box size
        x_min_box_B13= round( - 3* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B13= round( - 3* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B13= round( x_min_box_B13 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B13= round( y_min_box_B13 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const) 
        
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B13 = RectangleMesh(Point(x_min_box_B13, y_min_box_B13), Point(x_max_box_B13, y_max_box_B13), nx_B13, ny_B13)                                       
        V_B13 = FunctionSpace(mesh_B13, 'Lagrange', 1)        
        
        n_mesh_B13 = V_B13.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B13 = mesh_B13.geometry().dim()                                                        
        dof_coordinates_B13 = V_B13.tabulate_dof_coordinates()                   
        dof_coordinates_B13.resize((n_mesh_B13, d_mesh_B13))                                                   
        dof_x_B13 = dof_coordinates_B13[:, 0]                                                    
        dof_y_B13 = dof_coordinates_B13[:, 1]    
        
        
        
        
        ############################################################### 
        ## B1 mesh forming
        ############################################################### 
            
        ## mesh box size
        x_min_box_B1= round( - 7* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_min_box_B1= round( - 7* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B1= round( x_min_box_B1 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B1= round( y_min_box_B1 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const) 
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B1 = RectangleMesh(Point(x_min_box_B1, y_min_box_B1), Point(x_max_box_B1, y_max_box_B1), nx_B13, ny_B13)                                       
        V_B1 = FunctionSpace(mesh_B1, 'Lagrange', 1)        
        
        n_mesh_B1 = V_B1.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B1 = mesh_B1.geometry().dim()                                                        
        dof_coordinates_B1 = V_B1.tabulate_dof_coordinates()                   
        dof_coordinates_B1.resize((n_mesh_B1, d_mesh_B1))                                                   
        dof_x_B1 = dof_coordinates_B1[:, 0]                                                    
        dof_y_B1 = dof_coordinates_B1[:, 1]  
        
        
        
        ############################################################### 
        ## B2 mesh forming
        ############################################################### 
            
        ## mesh box size
        x_min_box_B2= round( - 5* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B2= round( - 7* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B2= round( x_min_box_B2 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B2= round( y_min_box_B2 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const) 
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B2 = RectangleMesh(Point(x_min_box_B2, y_min_box_B2), Point(x_max_box_B2, y_max_box_B2), nx_B13, ny_B13)                                       
        V_B2 = FunctionSpace(mesh_B2, 'Lagrange', 1)        
        
        n_mesh_B2 = V_B2.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B2 = mesh_B2.geometry().dim()                                                        
        dof_coordinates_B2 = V_B2.tabulate_dof_coordinates()                   
        dof_coordinates_B2.resize((n_mesh_B2, d_mesh_B2))                                                   
        dof_x_B2 = dof_coordinates_B2[:, 0]                                                    
        dof_y_B2 = dof_coordinates_B2[:, 1] 
        
        
        
        ############################################################### 
        ## B3 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B3= round( - 3* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B3= round( - 7* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B3= round( x_min_box_B3 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B3= round( y_min_box_B3 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B3 = RectangleMesh(Point(x_min_box_B3, y_min_box_B3), Point(x_max_box_B3, y_max_box_B3), nx_B13, ny_B13)                                       
        V_B3 = FunctionSpace(mesh_B3, 'Lagrange', 1)        
        
        n_mesh_B3 = V_B3.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B3 = mesh_B3.geometry().dim()                                                        
        dof_coordinates_B3 = V_B3.tabulate_dof_coordinates()                   
        dof_coordinates_B3.resize((n_mesh_B3, d_mesh_B3))                                                   
        dof_x_B3 = dof_coordinates_B3[:, 0]                                                    
        dof_y_B3 = dof_coordinates_B3[:, 1] 
        
        
        ############################################################### 
        ## B4 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B4= round( - 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B4= round( - 7* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B4= round( x_min_box_B4 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B4= round( y_min_box_B4 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B4 = RectangleMesh(Point(x_min_box_B4, y_min_box_B4), Point(x_max_box_B4, y_max_box_B4), nx_B13, ny_B13)                                       
        V_B4 = FunctionSpace(mesh_B4, 'Lagrange', 1)        
        
        n_mesh_B4 = V_B4.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B4 = mesh_B4.geometry().dim()                                                        
        dof_coordinates_B4 = V_B4.tabulate_dof_coordinates()                   
        dof_coordinates_B4.resize((n_mesh_B4, d_mesh_B4))                                                   
        dof_x_B4 = dof_coordinates_B4[:, 0]                                                    
        dof_y_B4 = dof_coordinates_B4[:, 1]  
        
        
        ############################################################### 
        ## B5 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B5= round( 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const) 
        y_min_box_B5= round( - 7* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B5= round( x_min_box_B5 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B5= round( y_min_box_B5 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B5 = RectangleMesh(Point(x_min_box_B5, y_min_box_B5), Point(x_max_box_B5, y_max_box_B5), nx_B13, ny_B13)                                       
        V_B5 = FunctionSpace(mesh_B5, 'Lagrange', 1)        
        
        n_mesh_B5 = V_B5.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B5 = mesh_B5.geometry().dim()                                                        
        dof_coordinates_B5 = V_B5.tabulate_dof_coordinates()                   
        dof_coordinates_B5.resize((n_mesh_B5, d_mesh_B5))                                                   
        dof_x_B5 = dof_coordinates_B5[:, 0]                                                    
        dof_y_B5 = dof_coordinates_B5[:, 1]  
        
        
        ############################################################### 
        ## B6 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B6= round( - 7* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B6= round( - 5* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B6= round( x_min_box_B6 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B6= round( y_min_box_B6 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B6 = RectangleMesh(Point(x_min_box_B6, y_min_box_B6), Point(x_max_box_B6, y_max_box_B6), nx_B13, ny_B13)                                       
        V_B6 = FunctionSpace(mesh_B6, 'Lagrange', 1)        
        
        n_mesh_B6 = V_B6.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B6 = mesh_B6.geometry().dim()                                                        
        dof_coordinates_B6 = V_B6.tabulate_dof_coordinates()                   
        dof_coordinates_B6.resize((n_mesh_B6, d_mesh_B6))                                                   
        dof_x_B6 = dof_coordinates_B6[:, 0]                                                    
        dof_y_B6 = dof_coordinates_B6[:, 1]  
        
        
        ############################################################### 
        ## B7 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B7= round( - 5* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B7= round( - 5* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const) 
        
        x_max_box_B7= round( x_min_box_B7 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B7= round( y_min_box_B7 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B7 = RectangleMesh(Point(x_min_box_B7, y_min_box_B7), Point(x_max_box_B7, y_max_box_B7), nx_B13, ny_B13)                                       
        V_B7 = FunctionSpace(mesh_B7, 'Lagrange', 1)        
        
        n_mesh_B7 = V_B7.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B7 = mesh_B7.geometry().dim()                                                        
        dof_coordinates_B7 = V_B7.tabulate_dof_coordinates()                   
        dof_coordinates_B7.resize((n_mesh_B7, d_mesh_B7))                                                   
        dof_x_B7 = dof_coordinates_B7[:, 0]                                                    
        dof_y_B7 = dof_coordinates_B7[:, 1]  
        
        
        ############################################################### 
        ## B8 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B8= round( - 3* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B8= round( - 5* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B8= round( x_min_box_B8 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B8= round( y_min_box_B8 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B8 = RectangleMesh(Point(x_min_box_B8, y_min_box_B8), Point(x_max_box_B8, y_max_box_B8), nx_B13, ny_B13)                                       
        V_B8 = FunctionSpace(mesh_B8, 'Lagrange', 1)        
        
        n_mesh_B8 = V_B8.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B8 = mesh_B8.geometry().dim()                                                        
        dof_coordinates_B8 = V_B8.tabulate_dof_coordinates()                   
        dof_coordinates_B8.resize((n_mesh_B8, d_mesh_B8))                                                   
        dof_x_B8 = dof_coordinates_B8[:, 0]                                                    
        dof_y_B8 = dof_coordinates_B8[:, 1]  
        
        
        ############################################################### 
        ## B9 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B9= round( - 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B9= round( - 5* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const) 
        
        x_max_box_B9= round( x_min_box_B9 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B9= round( y_min_box_B9 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B9 = RectangleMesh(Point(x_min_box_B9, y_min_box_B9), Point(x_max_box_B9, y_max_box_B9), nx_B13, ny_B13)                                       
        V_B9 = FunctionSpace(mesh_B9, 'Lagrange', 1)        
        
        n_mesh_B9 = V_B9.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B9 = mesh_B9.geometry().dim()                                                        
        dof_coordinates_B9 = V_B9.tabulate_dof_coordinates()                   
        dof_coordinates_B9.resize((n_mesh_B9, d_mesh_B9))                                                   
        dof_x_B9 = dof_coordinates_B9[:, 0]                                                    
        dof_y_B9 = dof_coordinates_B9[:, 1]  
        
        
        ############################################################### 
        ## B10 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B10= round( 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B10= round( - 5* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const) 
        
        x_max_box_B10= round( x_min_box_B10 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B10= round( y_min_box_B10 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B10 = RectangleMesh(Point(x_min_box_B10, y_min_box_B10), Point(x_max_box_B10, y_max_box_B10), nx_B13, ny_B13)                                       
        V_B10 = FunctionSpace(mesh_B10, 'Lagrange', 1)        
        
        n_mesh_B10 = V_B10.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B10 = mesh_B10.geometry().dim()                                                        
        dof_coordinates_B10 = V_B10.tabulate_dof_coordinates()                   
        dof_coordinates_B10.resize((n_mesh_B10, d_mesh_B10))                                                   
        dof_x_B10 = dof_coordinates_B10[:, 0]                                                    
        dof_y_B10 = dof_coordinates_B10[:, 1]  
        
        
        ############################################################### 
        ## B11 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B11= round( - 7* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B11= round( - 3* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B11= round( x_min_box_B11 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B11= round( y_min_box_B11 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B11 = RectangleMesh(Point(x_min_box_B11, y_min_box_B11), Point(x_max_box_B11, y_max_box_B11), nx_B13, ny_B13)                                       
        V_B11 = FunctionSpace(mesh_B11, 'Lagrange', 1)        
        
        n_mesh_B11 = V_B11.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B11 = mesh_B11.geometry().dim()                                                        
        dof_coordinates_B11 = V_B11.tabulate_dof_coordinates()                   
        dof_coordinates_B11.resize((n_mesh_B11, d_mesh_B11))                                                   
        dof_x_B11 = dof_coordinates_B11[:, 0]                                                    
        dof_y_B11 = dof_coordinates_B11[:, 1]  
        
        
        ############################################################### 
        ## B12 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B12= round( - 5* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B12= round( - 3* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B12= round( x_min_box_B12 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B12= round( y_min_box_B12 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B12 = RectangleMesh(Point(x_min_box_B12, y_min_box_B12), Point(x_max_box_B12, y_max_box_B12), nx_B13, ny_B13)                                       
        V_B12 = FunctionSpace(mesh_B12, 'Lagrange', 1)        
        
        n_mesh_B12 = V_B12.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B12 = mesh_B12.geometry().dim()                                                        
        dof_coordinates_B12 = V_B12.tabulate_dof_coordinates()                   
        dof_coordinates_B12.resize((n_mesh_B12, d_mesh_B12))                                                   
        dof_x_B12 = dof_coordinates_B12[:, 0]                                                    
        dof_y_B12 = dof_coordinates_B12[:, 1]  
         
        
        
        ############################################################### 
        ## B14 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B14= round( - 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B14= round( - 3* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B14= round( x_min_box_B14 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B14= round( y_min_box_B14 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B14 = RectangleMesh(Point(x_min_box_B14, y_min_box_B14), Point(x_max_box_B14, y_max_box_B14), nx_B13, ny_B13)                                       
        V_B14 = FunctionSpace(mesh_B14, 'Lagrange', 1)        
        
        n_mesh_B14 = V_B14.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B14 = mesh_B14.geometry().dim()                                                        
        dof_coordinates_B14 = V_B14.tabulate_dof_coordinates()                   
        dof_coordinates_B14.resize((n_mesh_B14, d_mesh_B14))                                                   
        dof_x_B14 = dof_coordinates_B14[:, 0]                                                    
        dof_y_B14 = dof_coordinates_B14[:, 1]  
        
        
        ############################################################### 
        ## B15 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B15= round(  1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B15= round( - 3* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B15= round( x_min_box_B15 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B15= round( y_min_box_B15 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B15 = RectangleMesh(Point(x_min_box_B15, y_min_box_B15), Point(x_max_box_B15, y_max_box_B15), nx_B13, ny_B13)                                       
        V_B15 = FunctionSpace(mesh_B15, 'Lagrange', 1)        
        
        n_mesh_B15 = V_B15.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B15 = mesh_B15.geometry().dim()                                                        
        dof_coordinates_B15 = V_B15.tabulate_dof_coordinates()                   
        dof_coordinates_B15.resize((n_mesh_B15, d_mesh_B15))                                                   
        dof_x_B15 = dof_coordinates_B15[:, 0]                                                    
        dof_y_B15 = dof_coordinates_B15[:, 1]  
        
        
        ############################################################### 
        ## B16 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B16= round( - 7* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B16= round( - 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B16= round( x_min_box_B16 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B16= round( y_min_box_B16 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B16 = RectangleMesh(Point(x_min_box_B16, y_min_box_B16), Point(x_max_box_B16, y_max_box_B16), nx_B13, ny_B13)                                       
        V_B16 = FunctionSpace(mesh_B16, 'Lagrange', 1)        
        
        n_mesh_B16 = V_B16.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B16 = mesh_B16.geometry().dim()                                                        
        dof_coordinates_B16 = V_B16.tabulate_dof_coordinates()                   
        dof_coordinates_B16.resize((n_mesh_B16, d_mesh_B16))                                                   
        dof_x_B16 = dof_coordinates_B16[:, 0]                                                    
        dof_y_B16 = dof_coordinates_B16[:, 1]  
        
        
        ############################################################### 
        ## B17 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B17= round( - 5* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B17= round( - 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B17= round( x_min_box_B17 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B17= round( y_min_box_B17 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B17 = RectangleMesh(Point(x_min_box_B17, y_min_box_B17), Point(x_max_box_B17, y_max_box_B17), nx_B13, ny_B13)                                       
        V_B17 = FunctionSpace(mesh_B17, 'Lagrange', 1)        
        
        n_mesh_B17 = V_B17.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B17 = mesh_B17.geometry().dim()                                                        
        dof_coordinates_B17 = V_B17.tabulate_dof_coordinates()                   
        dof_coordinates_B17.resize((n_mesh_B17, d_mesh_B17))                                                   
        dof_x_B17 = dof_coordinates_B17[:, 0]                                                    
        dof_y_B17 = dof_coordinates_B17[:, 1]  
        
        
        ############################################################### 
        ## B18 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B18= round( - 3* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B18= round( - 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B18= round( x_min_box_B18 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B18= round( y_min_box_B18 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B18 = RectangleMesh(Point(x_min_box_B18, y_min_box_B18), Point(x_max_box_B18, y_max_box_B18), nx_B13, ny_B13)                                       
        V_B18 = FunctionSpace(mesh_B18, 'Lagrange', 1)        
        
        n_mesh_B18 = V_B18.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B18 = mesh_B18.geometry().dim()                                                        
        dof_coordinates_B18 = V_B18.tabulate_dof_coordinates()                   
        dof_coordinates_B18.resize((n_mesh_B18, d_mesh_B18))                                                   
        dof_x_B18 = dof_coordinates_B18[:, 0]                                                    
        dof_y_B18 = dof_coordinates_B18[:, 1]  
        
        
        ############################################################### 
        ## B19 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B19= round( - 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B19= round( - 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B19= round( x_min_box_B19 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B19= round( y_min_box_B19 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) )  , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B19 = RectangleMesh(Point(x_min_box_B19, y_min_box_B19), Point(x_max_box_B19, y_max_box_B19), nx_B13, ny_B13)                                       
        V_B19 = FunctionSpace(mesh_B19, 'Lagrange', 1)        
        
        n_mesh_B19 = V_B19.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B19 = mesh_B19.geometry().dim()                                                        
        dof_coordinates_B19 = V_B19.tabulate_dof_coordinates()                   
        dof_coordinates_B19.resize((n_mesh_B19, d_mesh_B19))                                                   
        dof_x_B19 = dof_coordinates_B19[:, 0]                                                    
        dof_y_B19 = dof_coordinates_B19[:, 1]  
        
        
        ############################################################### 
        ## B20 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B20=  round( 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B20= round( - 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        
        x_max_box_B20= round( x_min_box_B20 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B20= round( y_min_box_B20 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B20 = RectangleMesh(Point(x_min_box_B20, y_min_box_B20), Point(x_max_box_B20, y_max_box_B20), nx_B13, ny_B13)                                       
        V_B20 = FunctionSpace(mesh_B20, 'Lagrange', 1)        
        
        n_mesh_B20 = V_B20.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B20 = mesh_B20.geometry().dim()                                                        
        dof_coordinates_B20 = V_B20.tabulate_dof_coordinates()                   
        dof_coordinates_B20.resize((n_mesh_B20, d_mesh_B20))                                                   
        dof_x_B20 = dof_coordinates_B20[:, 0]                                                    
        dof_y_B20 = dof_coordinates_B20[:, 1]  
        
        
        ############################################################### 
        ## B21 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B21= round( - 7* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B21=  round( 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B21= round( x_min_box_B21 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_max_box_B21= round( y_min_box_B21 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B21 = RectangleMesh(Point(x_min_box_B21, y_min_box_B21), Point(x_max_box_B21, y_max_box_B21), nx_B13, ny_B13)                                       
        V_B21 = FunctionSpace(mesh_B21, 'Lagrange', 1)        
        
        n_mesh_B21 = V_B21.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B21 = mesh_B21.geometry().dim()                                                        
        dof_coordinates_B21 = V_B21.tabulate_dof_coordinates()                   
        dof_coordinates_B21.resize((n_mesh_B21, d_mesh_B21))                                                   
        dof_x_B21 = dof_coordinates_B21[:, 0]                                                    
        dof_y_B21 = dof_coordinates_B21[:, 1]  
        
        
        ############################################################### 
        ## B22 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B22= round( - 5* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B22=  round( 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B22= round( x_min_box_B22 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B22= round( y_min_box_B22 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B22 = RectangleMesh(Point(x_min_box_B22, y_min_box_B22), Point(x_max_box_B22, y_max_box_B22), nx_B13, ny_B13)                                       
        V_B22 = FunctionSpace(mesh_B22, 'Lagrange', 1)        
        
        n_mesh_B22 = V_B22.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B22 = mesh_B22.geometry().dim()                                                        
        dof_coordinates_B22 = V_B22.tabulate_dof_coordinates()                   
        dof_coordinates_B22.resize((n_mesh_B22, d_mesh_B22))                                                   
        dof_x_B22 = dof_coordinates_B22[:, 0]                                                    
        dof_y_B22 = dof_coordinates_B22[:, 1]  
        
        
        ############################################################### 
        ## B23 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B23= round( - 3* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B23=  round( 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B23= round( x_min_box_B23 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B23= round( y_min_box_B23 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B23 = RectangleMesh(Point(x_min_box_B23, y_min_box_B23), Point(x_max_box_B23, y_max_box_B23), nx_B13, ny_B13)                                       
        V_B23 = FunctionSpace(mesh_B23, 'Lagrange', 1)        
        
        n_mesh_B23 = V_B23.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B23 = mesh_B23.geometry().dim()                                                        
        dof_coordinates_B23 = V_B23.tabulate_dof_coordinates()                   
        dof_coordinates_B23.resize((n_mesh_B23, d_mesh_B23))                                                   
        dof_x_B23 = dof_coordinates_B23[:, 0]                                                    
        dof_y_B23 = dof_coordinates_B23[:, 1]  
        
        
        ############################################################### 
        ## B24 mesh forming
        ############################################################### 
            
        ## mesh box size
        x_min_box_B24= round( - 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B24=  round( 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B24= round( x_min_box_B24 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B24= round( y_min_box_B24 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B24 = RectangleMesh(Point(x_min_box_B24, y_min_box_B24), Point(x_max_box_B24, y_max_box_B24), nx_B13, ny_B13)                                       
        V_B24 = FunctionSpace(mesh_B24, 'Lagrange', 1)        
        
        n_mesh_B24 = V_B24.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B24 = mesh_B24.geometry().dim()                                                        
        dof_coordinates_B24 = V_B24.tabulate_dof_coordinates()                   
        dof_coordinates_B24.resize((n_mesh_B24, d_mesh_B24))                                                   
        dof_x_B24 = dof_coordinates_B24[:, 0]                                                    
        dof_y_B24 = dof_coordinates_B24[:, 1]  
        
        
        ############################################################### 
        ## B25 mesh forming
        ############################################################### 
    
        ## mesh box size
        x_min_box_B25=  round( 1* (lambda_1_stretch[i1])  *  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        y_min_box_B25=  round( 1* (lambda_2_stretch[i2])*  ( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        x_max_box_B25= round( x_min_box_B25 +  6* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_B25= round( y_min_box_B25 +  6* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ) , round_const)
        
        
        ## Create mesh and define function space and dof coordinates
        mesh_B25 = RectangleMesh(Point(x_min_box_B25, y_min_box_B25), Point(x_max_box_B25, y_max_box_B25), nx_B13, ny_B13)                                       
        V_B25 = FunctionSpace(mesh_B25, 'Lagrange', 1)        
        
        n_mesh_B25 = V_B25.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B25 = mesh_B25.geometry().dim()                                                        
        dof_coordinates_B25 = V_B25.tabulate_dof_coordinates()                   
        dof_coordinates_B25.resize((n_mesh_B25, d_mesh_B25))                                                   
        dof_x_B25 = dof_coordinates_B25[:, 0]                                                    
        dof_y_B25 = dof_coordinates_B25[:, 1]  
        
        
        
       
        ############################################################### 
        ## Biggest mesh forming
        ############################################################### 
    
        # mesh size
        nx_big =    7* int( round( no_elements_per_lengthscale*(lambda_1_stretch[i1])   *2*  (  1/(2*sqrt(2))   ) ) )  #no of elements along x
        ny_big =    7* int( round( no_elements_per_lengthscale*(lambda_1_stretch[i1]) *2*  (  1/(2*sqrt(2))   ) ) )    #no of elements along y
        
        # nx_big= round( (7/3)*nx_V )
        # ny_big= round( (7/3)*ny_V )      
        
        ## mesh box size
        x_min_box_big =  round( - 7* (lambda_1_stretch[i1])    *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_min_box_big =  round( - 7* (lambda_2_stretch[i2])  *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
    
        x_max_box_big =  round( x_min_box_big + 14* (lambda_1_stretch[i1])   *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        y_max_box_big =  round( y_min_box_big + 14* (lambda_2_stretch[i2]) *( (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) ), round_const)
        

        
        
        ## Create mesh and define function space and dof coordinates
        mesh_big = RectangleMesh(Point(x_min_box_big, y_min_box_big), Point(x_max_box_big, y_max_box_big), nx_big, ny_big)                                       
        V_big = FunctionSpace(mesh_big, 'Lagrange', 1 )        
    
        n_mesh_big = V_big.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_big = mesh_big.geometry().dim()                                                        
        dof_coordinates_big = V_big.tabulate_dof_coordinates()                   
        dof_coordinates_big.resize((n_mesh_big, d_mesh_big))                                                   
        dof_x_big = dof_coordinates_big[:, 0]                                                    
        dof_y_big = dof_coordinates_big[:, 1]   
    
    

    
        ############################################################### 
        ## submeshes forming
        ############################################################### 
    
    
    
        ############################################################### box B13_submesh (central box)
    
        ## B13_submesh box size
        x_min_box_B13_submesh= x_min_box_B13
        y_min_box_B13_submesh= y_min_box_B13
        
        x_max_box_B13_submesh= x_max_box_B13
        y_max_box_B13_submesh= y_max_box_B13
    
    
        class domain_B13_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B13_submesh + tol and x[0] > x_min_box_B13_submesh - tol and  x[1]< y_max_box_B13_submesh + tol  and  x[1] > y_min_box_B13_submesh - tol
        
        
        
        subdomain = MeshFunction("size_t", mesh_big, mesh_big.topology().dim(),0)
        # subdomain.set_all(0)
        
        domain_B13_submesh_temp = domain_B13_submesh()
        domain_B13_submesh_temp.mark(subdomain, 13 )
        
        mesh_B13_submesh = SubMesh(mesh_big, subdomain, 13)
        V_B13_submesh = FunctionSpace(mesh_B13_submesh,'Lagrange',1)
        # V_B13 = FunctionSpace(mesh_B13, 'Lagrange', 1, constrained_domain=PeriodicBoundary_B13())        
    
    
        n_mesh_B13_submesh = V_B13_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B13_submesh = mesh_B13_submesh.geometry().dim()                                                        
        dof_coordinates_B13_submesh = V_B13_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B13_submesh.resize((n_mesh_B13_submesh, d_mesh_B13_submesh))                                                   
        dof_x_B13_submesh = dof_coordinates_B13_submesh[:, 0]                                                    
        dof_y_B13_submesh = dof_coordinates_B13_submesh[:, 1]
    
    
    
        ############################################################### box B1_submesh 
    
        ## B1_submesh box size
        x_min_box_B1_submesh=  x_min_box_B1
        y_min_box_B1_submesh=  y_min_box_B1
        
        x_max_box_B1_submesh=  x_max_box_B1
        y_max_box_B1_submesh=  y_max_box_B1
    
    
        class domain_B1_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B1_submesh + tol and x[0] > x_min_box_B1_submesh - tol and  x[1]< y_max_box_B1_submesh + tol  and  x[1] > y_min_box_B1_submesh - tol
    
        

        domain_B1_submesh_temp = domain_B1_submesh()
        domain_B1_submesh_temp.mark(subdomain, 1 )
        
        mesh_B1_submesh = SubMesh(mesh_big, subdomain, 1)
        V_B1_submesh = FunctionSpace(mesh_B1_submesh,'Lagrange',1)
    
    
        n_mesh_B1_submesh = V_B1_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B1_submesh = mesh_B1_submesh.geometry().dim()                                                        
        dof_coordinates_B1_submesh = V_B1_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B1_submesh.resize((n_mesh_B1_submesh, d_mesh_B1_submesh))                                                   
        dof_x_B1_submesh = dof_coordinates_B1_submesh[:, 0]                                                    
        dof_y_B1_submesh = dof_coordinates_B1_submesh[:, 1]
    
    
        # ############################################################### box B2_submesh 
    
        ## B2_submesh box size
        x_min_box_B2_submesh = x_min_box_B2
        y_min_box_B2_submesh = y_min_box_B2
        
        x_max_box_B2_submesh = x_max_box_B2
        y_max_box_B2_submesh = y_max_box_B2
    
    
        class domain_B2_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B2_submesh + tol and x[0] > x_min_box_B2_submesh - tol and  x[1]< y_max_box_B2_submesh + tol  and  x[1] > y_min_box_B2_submesh - tol
        
      
        
        domain_B2_submesh_temp = domain_B2_submesh()
        domain_B2_submesh_temp.mark(subdomain, 2 )
        
        mesh_B2_submesh = SubMesh(mesh_big, subdomain, 2)
        V_B2_submesh = FunctionSpace(mesh_B2_submesh,'Lagrange',1)
    
    
        n_mesh_B2_submesh = V_B2_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B2_submesh = mesh_B2_submesh.geometry().dim()                                                        
        dof_coordinates_B2_submesh = V_B2_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B2_submesh.resize((n_mesh_B2_submesh, d_mesh_B2_submesh))                                                   
        dof_x_B2_submesh = dof_coordinates_B2_submesh[:, 0]                                                    
        dof_y_B2_submesh = dof_coordinates_B2_submesh[:, 1]
    
    
        ############################################################### box B3_submesh
    
        ## B3_submesh box size
        x_min_box_B3_submesh= x_min_box_B3
        y_min_box_B3_submesh= y_min_box_B3
        
        x_max_box_B3_submesh= x_max_box_B3
        y_max_box_B3_submesh= y_max_box_B3
    
    
        class domain_B3_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B3_submesh + tol and x[0] > x_min_box_B3_submesh - tol and  x[1]< y_max_box_B3_submesh + tol  and  x[1] > y_min_box_B3_submesh - tol
        
         
        
        domain_B3_submesh_temp = domain_B3_submesh()
        domain_B3_submesh_temp.mark(subdomain, 3 )
        
        mesh_B3_submesh = SubMesh(mesh_big, subdomain, 3)
        V_B3_submesh = FunctionSpace(mesh_B3_submesh,'Lagrange',1)
    
    
        n_mesh_B3_submesh = V_B3_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B3_submesh = mesh_B3_submesh.geometry().dim()                                                        
        dof_coordinates_B3_submesh = V_B3_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B3_submesh.resize((n_mesh_B3_submesh, d_mesh_B3_submesh))                                                   
        dof_x_B3_submesh = dof_coordinates_B3_submesh[:, 0]                                                    
        dof_y_B3_submesh = dof_coordinates_B3_submesh[:, 1]
    
    
        ############################################################### box B4_submesh 
    
        ## B4_submesh box size
        x_min_box_B4_submesh=  x_min_box_B4
        y_min_box_B4_submesh=  y_min_box_B4
        
        x_max_box_B4_submesh=  x_max_box_B4
        y_max_box_B4_submesh=  y_max_box_B4
    
    
        class domain_B4_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B4_submesh + tol and x[0] > x_min_box_B4_submesh - tol and  x[1]< y_max_box_B4_submesh + tol  and  x[1] > y_min_box_B4_submesh - tol
    
        
        domain_B4_submesh_temp = domain_B4_submesh()
        domain_B4_submesh_temp.mark(subdomain, 4 )
        
        mesh_B4_submesh = SubMesh(mesh_big, subdomain, 4)
        V_B4_submesh = FunctionSpace(mesh_B4_submesh,'Lagrange',1)
    
    
        n_mesh_B4_submesh = V_B4_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B4_submesh = mesh_B4_submesh.geometry().dim()                                                        
        dof_coordinates_B4_submesh = V_B4_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B4_submesh.resize((n_mesh_B4_submesh, d_mesh_B4_submesh))                                                   
        dof_x_B4_submesh = dof_coordinates_B4_submesh[:, 0]                                                    
        dof_y_B4_submesh = dof_coordinates_B4_submesh[:, 1]
    
    
    
        ############################################################### box B5_submesh 
    
        ## B5_submesh box size
        x_min_box_B5_submesh=  x_min_box_B5
        y_min_box_B5_submesh=  y_min_box_B5
        
        x_max_box_B5_submesh=  x_max_box_B5
        y_max_box_B5_submesh=  y_max_box_B5
    
    
        class domain_B5_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B5_submesh + tol and x[0] > x_min_box_B5_submesh - tol and  x[1]< y_max_box_B5_submesh + tol  and  x[1] > y_min_box_B5_submesh - tol
    
        
        domain_B5_submesh_temp = domain_B5_submesh()
        domain_B5_submesh_temp.mark(subdomain, 5 )
        
        mesh_B5_submesh = SubMesh(mesh_big, subdomain, 5)
        V_B5_submesh = FunctionSpace(mesh_B5_submesh,'Lagrange',1)
    
    
        n_mesh_B5_submesh = V_B5_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B5_submesh = mesh_B5_submesh.geometry().dim()                                                        
        dof_coordinates_B5_submesh = V_B5_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B5_submesh.resize((n_mesh_B5_submesh, d_mesh_B5_submesh))                                                   
        dof_x_B5_submesh = dof_coordinates_B5_submesh[:, 0]                                                    
        dof_y_B5_submesh = dof_coordinates_B5_submesh[:, 1]
        
        
        ############################################################### box B6_submesh 
    
        ## B6_submesh box size
        x_min_box_B6_submesh=  x_min_box_B6
        y_min_box_B6_submesh=  y_min_box_B6
        
        x_max_box_B6_submesh=  x_max_box_B6
        y_max_box_B6_submesh=  y_max_box_B6
    
    
        class domain_B6_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B6_submesh + tol and x[0] > x_min_box_B6_submesh - tol and  x[1]< y_max_box_B6_submesh + tol  and  x[1] > y_min_box_B6_submesh - tol
    
        
        domain_B6_submesh_temp = domain_B6_submesh()
        domain_B6_submesh_temp.mark(subdomain, 6 )
        
        mesh_B6_submesh = SubMesh(mesh_big, subdomain, 6)
        V_B6_submesh = FunctionSpace(mesh_B6_submesh,'Lagrange',1)
    
    
        n_mesh_B6_submesh = V_B6_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B6_submesh = mesh_B6_submesh.geometry().dim()                                                        
        dof_coordinates_B6_submesh = V_B6_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B6_submesh.resize((n_mesh_B6_submesh, d_mesh_B6_submesh))                                                   
        dof_x_B6_submesh = dof_coordinates_B6_submesh[:, 0]                                                    
        dof_y_B6_submesh = dof_coordinates_B6_submesh[:, 1]
        
        
        ############################################################### box B7_submesh 
    
        ## B7_submesh box size
        x_min_box_B7_submesh=  x_min_box_B7
        y_min_box_B7_submesh=  y_min_box_B7
        
        x_max_box_B7_submesh=  x_max_box_B7
        y_max_box_B7_submesh=  y_max_box_B7
    
    
        class domain_B7_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B7_submesh + tol and x[0] > x_min_box_B7_submesh - tol and  x[1]< y_max_box_B7_submesh + tol  and  x[1] > y_min_box_B7_submesh - tol
    
        
        domain_B7_submesh_temp = domain_B7_submesh()
        domain_B7_submesh_temp.mark(subdomain, 7 )
        
        mesh_B7_submesh = SubMesh(mesh_big, subdomain, 7)
        V_B7_submesh = FunctionSpace(mesh_B7_submesh,'Lagrange',1)
    
    
        n_mesh_B7_submesh = V_B7_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B7_submesh = mesh_B7_submesh.geometry().dim()                                                        
        dof_coordinates_B7_submesh = V_B7_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B7_submesh.resize((n_mesh_B7_submesh, d_mesh_B7_submesh))                                                   
        dof_x_B7_submesh = dof_coordinates_B7_submesh[:, 0]                                                    
        dof_y_B7_submesh = dof_coordinates_B7_submesh[:, 1]
        
        
        ############################################################### box B8_submesh 
    
        ## B8_submesh box size
        x_min_box_B8_submesh=  x_min_box_B8
        y_min_box_B8_submesh=  y_min_box_B8
        
        x_max_box_B8_submesh=  x_max_box_B8
        y_max_box_B8_submesh=  y_max_box_B8
    
    
        class domain_B8_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B8_submesh + tol and x[0] > x_min_box_B8_submesh - tol and  x[1]< y_max_box_B8_submesh + tol  and  x[1] > y_min_box_B8_submesh - tol
    
        
        domain_B8_submesh_temp = domain_B8_submesh()
        domain_B8_submesh_temp.mark(subdomain, 8 )
        
        mesh_B8_submesh = SubMesh(mesh_big, subdomain, 8)
        V_B8_submesh = FunctionSpace(mesh_B8_submesh,'Lagrange',1)
    
    
        n_mesh_B8_submesh = V_B8_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B8_submesh = mesh_B8_submesh.geometry().dim()                                                        
        dof_coordinates_B8_submesh = V_B8_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B8_submesh.resize((n_mesh_B8_submesh, d_mesh_B8_submesh))                                                   
        dof_x_B8_submesh = dof_coordinates_B8_submesh[:, 0]                                                    
        dof_y_B8_submesh = dof_coordinates_B8_submesh[:, 1]
        
        
        ############################################################### box B9_submesh 
    
        ## B9_submesh box size
        x_min_box_B9_submesh=  x_min_box_B9
        y_min_box_B9_submesh=  y_min_box_B9
        
        x_max_box_B9_submesh=  x_max_box_B9
        y_max_box_B9_submesh=  y_max_box_B9
    
    
        class domain_B9_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B9_submesh + tol and x[0] > x_min_box_B9_submesh - tol and  x[1]< y_max_box_B9_submesh + tol  and  x[1] > y_min_box_B9_submesh - tol
    
        
        domain_B9_submesh_temp = domain_B9_submesh()
        domain_B9_submesh_temp.mark(subdomain, 9 )
        
        mesh_B9_submesh = SubMesh(mesh_big, subdomain, 9)
        V_B9_submesh = FunctionSpace(mesh_B9_submesh,'Lagrange',1)
    
    
        n_mesh_B9_submesh = V_B9_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B9_submesh = mesh_B9_submesh.geometry().dim()                                                        
        dof_coordinates_B9_submesh = V_B9_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B9_submesh.resize((n_mesh_B9_submesh, d_mesh_B9_submesh))                                                   
        dof_x_B9_submesh = dof_coordinates_B9_submesh[:, 0]                                                    
        dof_y_B9_submesh = dof_coordinates_B9_submesh[:, 1]
        
        
        ############################################################### box B10_submesh 
    
        ## B10_submesh box size
        x_min_box_B10_submesh=  x_min_box_B10
        y_min_box_B10_submesh=  y_min_box_B10
        
        x_max_box_B10_submesh=  x_max_box_B10
        y_max_box_B10_submesh=  y_max_box_B10
    
    
        class domain_B10_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B10_submesh + tol and x[0] > x_min_box_B10_submesh - tol and  x[1]< y_max_box_B10_submesh + tol  and  x[1] > y_min_box_B10_submesh - tol
        
        domain_B10_submesh_temp = domain_B10_submesh()
        domain_B10_submesh_temp.mark(subdomain, 10 )
        
        mesh_B10_submesh = SubMesh(mesh_big, subdomain, 10)
        V_B10_submesh = FunctionSpace(mesh_B10_submesh,'Lagrange',1)
    
    
        n_mesh_B10_submesh = V_B10_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B10_submesh = mesh_B10_submesh.geometry().dim()                                                        
        dof_coordinates_B10_submesh = V_B10_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B10_submesh.resize((n_mesh_B10_submesh, d_mesh_B10_submesh))                                                   
        dof_x_B10_submesh = dof_coordinates_B10_submesh[:, 0]                                                    
        dof_y_B10_submesh = dof_coordinates_B10_submesh[:, 1]
        
        
        ############################################################### box B11_submesh 
    
        ## B11_submesh box size
        x_min_box_B11_submesh=  x_min_box_B11
        y_min_box_B11_submesh=  y_min_box_B11
        
        x_max_box_B11_submesh=  x_max_box_B11
        y_max_box_B11_submesh=  y_max_box_B11
    
    
        class domain_B11_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B11_submesh + tol and x[0] > x_min_box_B11_submesh - tol and  x[1]< y_max_box_B11_submesh + tol  and  x[1] > y_min_box_B11_submesh - tol
    
        
        domain_B11_submesh_temp = domain_B11_submesh()
        domain_B11_submesh_temp.mark(subdomain, 11 )
        
        mesh_B11_submesh = SubMesh(mesh_big, subdomain, 11)
        V_B11_submesh = FunctionSpace(mesh_B11_submesh,'Lagrange',1)
    
    
        n_mesh_B11_submesh = V_B11_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B11_submesh = mesh_B11_submesh.geometry().dim()                                                        
        dof_coordinates_B11_submesh = V_B11_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B11_submesh.resize((n_mesh_B11_submesh, d_mesh_B11_submesh))                                                   
        dof_x_B11_submesh = dof_coordinates_B11_submesh[:, 0]                                                    
        dof_y_B11_submesh = dof_coordinates_B11_submesh[:, 1]
        
        
        ############################################################### box B12_submesh 
    
        ## B12_submesh box size
        x_min_box_B12_submesh=  x_min_box_B12
        y_min_box_B12_submesh=  y_min_box_B12
        
        x_max_box_B12_submesh=  x_max_box_B12
        y_max_box_B12_submesh=  y_max_box_B12
    
    
        class domain_B12_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B12_submesh + tol and x[0] > x_min_box_B12_submesh - tol and  x[1]< y_max_box_B12_submesh + tol  and  x[1] > y_min_box_B12_submesh - tol
    
        
        domain_B12_submesh_temp = domain_B12_submesh()
        domain_B12_submesh_temp.mark(subdomain, 12 )
        
        mesh_B12_submesh = SubMesh(mesh_big, subdomain, 12)
        V_B12_submesh = FunctionSpace(mesh_B12_submesh,'Lagrange',1)
    
    
        n_mesh_B12_submesh = V_B12_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B12_submesh = mesh_B12_submesh.geometry().dim()                                                        
        dof_coordinates_B12_submesh = V_B12_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B12_submesh.resize((n_mesh_B12_submesh, d_mesh_B12_submesh))                                                   
        dof_x_B12_submesh = dof_coordinates_B12_submesh[:, 0]                                                    
        dof_y_B12_submesh = dof_coordinates_B12_submesh[:, 1]
        
        
        ############################################################### box B14_submesh 
    
        ## B14_submesh box size
        x_min_box_B14_submesh=  x_min_box_B14
        y_min_box_B14_submesh=  y_min_box_B14
        
        x_max_box_B14_submesh=  x_max_box_B14
        y_max_box_B14_submesh=  y_max_box_B14
    
    
        class domain_B14_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B14_submesh + tol and x[0] > x_min_box_B14_submesh - tol and  x[1]< y_max_box_B14_submesh + tol  and  x[1] > y_min_box_B14_submesh - tol
    
        
        domain_B14_submesh_temp = domain_B14_submesh()
        domain_B14_submesh_temp.mark(subdomain, 14 )
        
        mesh_B14_submesh = SubMesh(mesh_big, subdomain, 14)
        V_B14_submesh = FunctionSpace(mesh_B14_submesh,'Lagrange',1)
    
    
        n_mesh_B14_submesh = V_B14_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B14_submesh = mesh_B14_submesh.geometry().dim()                                                        
        dof_coordinates_B14_submesh = V_B14_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B14_submesh.resize((n_mesh_B14_submesh, d_mesh_B14_submesh))                                                   
        dof_x_B14_submesh = dof_coordinates_B14_submesh[:, 0]                                                    
        dof_y_B14_submesh = dof_coordinates_B14_submesh[:, 1]
        
        
        ############################################################### box B15_submesh 
    
        ## B15_submesh box size
        x_min_box_B15_submesh=  x_min_box_B15
        y_min_box_B15_submesh=  y_min_box_B15
        
        x_max_box_B15_submesh=  x_max_box_B15
        y_max_box_B15_submesh=  y_max_box_B15
    
    
        class domain_B15_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B15_submesh + tol and x[0] > x_min_box_B15_submesh - tol and  x[1]< y_max_box_B15_submesh + tol  and  x[1] > y_min_box_B15_submesh - tol
    
        
        domain_B15_submesh_temp = domain_B15_submesh()
        domain_B15_submesh_temp.mark(subdomain, 15 )
        
        mesh_B15_submesh = SubMesh(mesh_big, subdomain, 15)
        V_B15_submesh = FunctionSpace(mesh_B15_submesh,'Lagrange',1)
    
    
        n_mesh_B15_submesh = V_B15_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B15_submesh = mesh_B15_submesh.geometry().dim()                                                        
        dof_coordinates_B15_submesh = V_B15_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B15_submesh.resize((n_mesh_B15_submesh, d_mesh_B15_submesh))                                                   
        dof_x_B15_submesh = dof_coordinates_B15_submesh[:, 0]                                                    
        dof_y_B15_submesh = dof_coordinates_B15_submesh[:, 1]
        
        
        ############################################################### box B16_submesh 
    
        ## B16_submesh box size
        x_min_box_B16_submesh=  x_min_box_B16
        y_min_box_B16_submesh=  y_min_box_B16
        
        x_max_box_B16_submesh=  x_max_box_B16
        y_max_box_B16_submesh=  y_max_box_B16
    
    
        class domain_B16_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B16_submesh + tol and x[0] > x_min_box_B16_submesh - tol and  x[1]< y_max_box_B16_submesh + tol  and  x[1] > y_min_box_B16_submesh - tol
    
        
        domain_B16_submesh_temp = domain_B16_submesh()
        domain_B16_submesh_temp.mark(subdomain, 16 )
        
        mesh_B16_submesh = SubMesh(mesh_big, subdomain, 16)
        V_B16_submesh = FunctionSpace(mesh_B16_submesh,'Lagrange',1)
    
    
        n_mesh_B16_submesh = V_B16_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B16_submesh = mesh_B16_submesh.geometry().dim()                                                        
        dof_coordinates_B16_submesh = V_B16_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B16_submesh.resize((n_mesh_B16_submesh, d_mesh_B16_submesh))                                                   
        dof_x_B16_submesh = dof_coordinates_B16_submesh[:, 0]                                                    
        dof_y_B16_submesh = dof_coordinates_B16_submesh[:, 1]
        
        
        ############################################################### box B17_submesh 
    
        ## B17_submesh box size
        x_min_box_B17_submesh=  x_min_box_B17
        y_min_box_B17_submesh=  y_min_box_B17
        
        x_max_box_B17_submesh=  x_max_box_B17
        y_max_box_B17_submesh=  y_max_box_B17
    
    
        class domain_B17_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B17_submesh + tol and x[0] > x_min_box_B17_submesh - tol and  x[1]< y_max_box_B17_submesh + tol  and  x[1] > y_min_box_B17_submesh - tol
    
        
        domain_B17_submesh_temp = domain_B17_submesh()
        domain_B17_submesh_temp.mark(subdomain, 17 )
        
        mesh_B17_submesh = SubMesh(mesh_big, subdomain, 17)
        V_B17_submesh = FunctionSpace(mesh_B17_submesh,'Lagrange',1)
    
    
        n_mesh_B17_submesh = V_B17_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B17_submesh = mesh_B17_submesh.geometry().dim()                                                        
        dof_coordinates_B17_submesh = V_B17_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B17_submesh.resize((n_mesh_B17_submesh, d_mesh_B17_submesh))                                                   
        dof_x_B17_submesh = dof_coordinates_B17_submesh[:, 0]                                                    
        dof_y_B17_submesh = dof_coordinates_B17_submesh[:, 1]
        
        
        ############################################################### box B18_submesh 
    
        ## B18_submesh box size
        x_min_box_B18_submesh=  x_min_box_B18
        y_min_box_B18_submesh=  y_min_box_B18
        
        x_max_box_B18_submesh=  x_max_box_B18
        y_max_box_B18_submesh=  y_max_box_B18
    
    
        class domain_B18_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B18_submesh + tol and x[0] > x_min_box_B18_submesh - tol and  x[1]< y_max_box_B18_submesh + tol  and  x[1] > y_min_box_B18_submesh - tol
    
        
        domain_B18_submesh_temp = domain_B18_submesh()
        domain_B18_submesh_temp.mark(subdomain, 18 )
        
        mesh_B18_submesh = SubMesh(mesh_big, subdomain, 18)
        V_B18_submesh = FunctionSpace(mesh_B18_submesh,'Lagrange',1)
    
    
        n_mesh_B18_submesh = V_B18_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B18_submesh = mesh_B18_submesh.geometry().dim()                                                        
        dof_coordinates_B18_submesh = V_B18_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B18_submesh.resize((n_mesh_B18_submesh, d_mesh_B18_submesh))                                                   
        dof_x_B18_submesh = dof_coordinates_B18_submesh[:, 0]                                                    
        dof_y_B18_submesh = dof_coordinates_B18_submesh[:, 1]
        
        
        ############################################################### box B19_submesh 
    
        ## B19_submesh box size
        x_min_box_B19_submesh=  x_min_box_B19
        y_min_box_B19_submesh=  y_min_box_B19
        
        x_max_box_B19_submesh=  x_max_box_B19
        y_max_box_B19_submesh=  y_max_box_B19
    
    
        class domain_B19_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B19_submesh + tol and x[0] > x_min_box_B19_submesh - tol and  x[1]< y_max_box_B19_submesh + tol  and  x[1] > y_min_box_B19_submesh - tol
    
        
        domain_B19_submesh_temp = domain_B19_submesh()
        domain_B19_submesh_temp.mark(subdomain, 19 )
        
        mesh_B19_submesh = SubMesh(mesh_big, subdomain, 19)
        V_B19_submesh = FunctionSpace(mesh_B19_submesh,'Lagrange',1)
    
    
        n_mesh_B19_submesh = V_B19_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B19_submesh = mesh_B19_submesh.geometry().dim()                                                        
        dof_coordinates_B19_submesh = V_B19_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B19_submesh.resize((n_mesh_B19_submesh, d_mesh_B19_submesh))                                                   
        dof_x_B19_submesh = dof_coordinates_B19_submesh[:, 0]                                                    
        dof_y_B19_submesh = dof_coordinates_B19_submesh[:, 1]
        
        
        ############################################################### box B20_submesh 
    
        ## B20_submesh box size
        x_min_box_B20_submesh=  x_min_box_B20
        y_min_box_B20_submesh=  y_min_box_B20
        
        x_max_box_B20_submesh=  x_max_box_B20
        y_max_box_B20_submesh=  y_max_box_B20
    
    
        class domain_B20_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B20_submesh + tol and x[0] > x_min_box_B20_submesh - tol and  x[1]< y_max_box_B20_submesh + tol  and  x[1] > y_min_box_B20_submesh - tol
    
        
        domain_B20_submesh_temp = domain_B20_submesh()
        domain_B20_submesh_temp.mark(subdomain, 20 )
        
        mesh_B20_submesh = SubMesh(mesh_big, subdomain, 20)
        V_B20_submesh = FunctionSpace(mesh_B20_submesh,'Lagrange',1)
    
    
        n_mesh_B20_submesh = V_B20_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B20_submesh = mesh_B20_submesh.geometry().dim()                                                        
        dof_coordinates_B20_submesh = V_B20_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B20_submesh.resize((n_mesh_B20_submesh, d_mesh_B20_submesh))                                                   
        dof_x_B20_submesh = dof_coordinates_B20_submesh[:, 0]                                                    
        dof_y_B20_submesh = dof_coordinates_B20_submesh[:, 1]
        
        
        ############################################################### box B21_submesh 
    
        ## B21_submesh box size
        x_min_box_B21_submesh=  x_min_box_B21
        y_min_box_B21_submesh=  y_min_box_B21
        
        x_max_box_B21_submesh=  x_max_box_B21
        y_max_box_B21_submesh=  y_max_box_B21
    
    
        class domain_B21_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B21_submesh + tol and x[0] > x_min_box_B21_submesh - tol and  x[1]< y_max_box_B21_submesh + tol  and  x[1] > y_min_box_B21_submesh - tol
    
        
        domain_B21_submesh_temp = domain_B21_submesh()
        domain_B21_submesh_temp.mark(subdomain, 21 )
        
        mesh_B21_submesh = SubMesh(mesh_big, subdomain, 21)
        V_B21_submesh = FunctionSpace(mesh_B21_submesh,'Lagrange',1)
    
    
        n_mesh_B21_submesh = V_B21_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B21_submesh = mesh_B21_submesh.geometry().dim()                                                        
        dof_coordinates_B21_submesh = V_B21_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B21_submesh.resize((n_mesh_B21_submesh, d_mesh_B21_submesh))                                                   
        dof_x_B21_submesh = dof_coordinates_B21_submesh[:, 0]                                                    
        dof_y_B21_submesh = dof_coordinates_B21_submesh[:, 1]
        
        
        ############################################################### box B22_submesh 
    
        ## B22_submesh box size
        x_min_box_B22_submesh=  x_min_box_B22
        y_min_box_B22_submesh=  y_min_box_B22
        
        x_max_box_B22_submesh=  x_max_box_B22
        y_max_box_B22_submesh=  y_max_box_B22
    
    
        class domain_B22_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B22_submesh + tol and x[0] > x_min_box_B22_submesh - tol and  x[1]< y_max_box_B22_submesh + tol  and  x[1] > y_min_box_B22_submesh - tol
    
        
        domain_B22_submesh_temp = domain_B22_submesh()
        domain_B22_submesh_temp.mark(subdomain, 22 )
        
        mesh_B22_submesh = SubMesh(mesh_big, subdomain, 22)
        V_B22_submesh = FunctionSpace(mesh_B22_submesh,'Lagrange',1)
    
        n_mesh_B22_submesh = V_B22_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B22_submesh = mesh_B22_submesh.geometry().dim()                                                        
        dof_coordinates_B22_submesh = V_B22_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B22_submesh.resize((n_mesh_B22_submesh, d_mesh_B22_submesh))                                                   
        dof_x_B22_submesh = dof_coordinates_B22_submesh[:, 0]                                                    
        dof_y_B22_submesh = dof_coordinates_B22_submesh[:, 1]
        
        
        ############################################################### box B23_submesh 
    
        ## B23_submesh box size
        x_min_box_B23_submesh=  x_min_box_B23
        y_min_box_B23_submesh=  y_min_box_B23
        
        x_max_box_B23_submesh=  x_max_box_B23
        y_max_box_B23_submesh=  y_max_box_B23
    
    
        class domain_B23_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B23_submesh + tol and x[0] > x_min_box_B23_submesh - tol and  x[1]< y_max_box_B23_submesh + tol  and  x[1] > y_min_box_B23_submesh - tol
    
        
        domain_B23_submesh_temp = domain_B23_submesh()
        domain_B23_submesh_temp.mark(subdomain, 23 )
        
        mesh_B23_submesh = SubMesh(mesh_big, subdomain, 23)
        V_B23_submesh = FunctionSpace(mesh_B23_submesh,'Lagrange', 1)    
    
        n_mesh_B23_submesh = V_B23_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B23_submesh = mesh_B23_submesh.geometry().dim()                                                        
        dof_coordinates_B23_submesh = V_B23_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B23_submesh.resize((n_mesh_B23_submesh, d_mesh_B23_submesh))                                                   
        dof_x_B23_submesh = dof_coordinates_B23_submesh[:, 0]                                                    
        dof_y_B23_submesh = dof_coordinates_B23_submesh[:, 1]
        
        
        ############################################################### box B24_submesh 
    
        ## B24_submesh box size
        x_min_box_B24_submesh=  x_min_box_B24
        y_min_box_B24_submesh=  y_min_box_B24
        
        x_max_box_B24_submesh=  x_max_box_B24
        y_max_box_B24_submesh=  y_max_box_B24
    
    
        class domain_B24_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B24_submesh + tol and x[0] > x_min_box_B24_submesh - tol and  x[1]< y_max_box_B24_submesh + tol  and  x[1] > y_min_box_B24_submesh - tol
    
        
        domain_B24_submesh_temp = domain_B24_submesh()
        domain_B24_submesh_temp.mark(subdomain, 24 )
        
        mesh_B24_submesh = SubMesh(mesh_big, subdomain, 24)
        V_B24_submesh = FunctionSpace(mesh_B24_submesh,'Lagrange',1)    
    
        n_mesh_B24_submesh = V_B24_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B24_submesh = mesh_B24_submesh.geometry().dim()                                                        
        dof_coordinates_B24_submesh = V_B24_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B24_submesh.resize((n_mesh_B24_submesh, d_mesh_B24_submesh))                                                   
        dof_x_B24_submesh = dof_coordinates_B24_submesh[:, 0]                                                    
        dof_y_B24_submesh = dof_coordinates_B24_submesh[:, 1]
        
        
        ############################################################### box B25_submesh 
    
        ## B25_submesh box size
        x_min_box_B25_submesh=  x_min_box_B25
        y_min_box_B25_submesh=  y_min_box_B25
        
        x_max_box_B25_submesh=  x_max_box_B25
        y_max_box_B25_submesh=  y_max_box_B25
    
    
        class domain_B25_submesh(SubDomain):
            def inside(self,x,on_boundary):
                return x[0]< x_max_box_B25_submesh + tol and x[0] > x_min_box_B25_submesh - tol and  x[1]< y_max_box_B25_submesh + tol  and  x[1] > y_min_box_B25_submesh - tol
    
        
        domain_B25_submesh_temp = domain_B25_submesh()
        domain_B25_submesh_temp.mark(subdomain, 25 )
        
        mesh_B25_submesh = SubMesh(mesh_big, subdomain, 25)
        V_B25_submesh = FunctionSpace(mesh_B25_submesh,'Lagrange',1)    
    
        n_mesh_B25_submesh = V_B25_submesh.dim()   #no of dof points, n_mesh=(nx+1)*(ny+1)                                                                
        d_mesh_B25_submesh = mesh_B25_submesh.geometry().dim()                                                        
        dof_coordinates_B25_submesh = V_B25_submesh.tabulate_dof_coordinates()                   
        dof_coordinates_B25_submesh.resize((n_mesh_B25_submesh, d_mesh_B25_submesh))                                                   
        dof_x_B25_submesh = dof_coordinates_B25_submesh[:, 0]                                                    
        dof_y_B25_submesh = dof_coordinates_B25_submesh[:, 1]
    
    
        
    
    
    
        print('stage 0 complete')
    
    
    
    
    
    
        ################################################################################
        ## Function for computing q 
        
        def q_computation(X, w):
        
        
            #initialize q_numr and q_star_numr 
            q_numr=np.zeros((n_t, n_mesh)) 
            
            ##initial q and q_star at t=0
            q_n= Function(V) 
           
            
            ##initial condition for q
            q_0= Function(V)                                      
            q_0_numr =(a_chain*sqrt(N_chain))*(1/(sqrt(2*pi)*c_dirac))*np.exp( (-1/(2*c_dirac**2))* (   (dof_x- X[0])**2 + (dof_y- X[1])**2  ) )  
            q_0.vector()[:] = q_0_numr 
            
            ##initialize u value at t=0
            q_n.assign(q_0)
            
            #save q_0 in file
            vtkfile_q << (q_0, 0)
            
            #save q_0 in matrix
            q_numr[0][:]=q_0.vector()[:]
        
                
        
            ######## time stepping for computing q
            for n in range(1,n_t): 
                
                # print(n)
                
                t=dt*n #contor parameter value (s)
        
                ######################################### computing q

                
                #defining q, v
                q = TrialFunction(V)
                v = TestFunction(V)
                
                #a, L in fem weak form for fenics (with Crank-Nicolson time stepping)
                a= G_chain*(dt/2)*dot(grad(q), grad(v))*dx + q*v*dx  + dt*w*q*v*dx 
                L= ( q_n*v*dx - G_chain*(dt/2)*dot(grad(q_n), grad(v))*dx  - dt*w*q_n*v*dx)
                
                #solve variational problem
                q = Function(V)
                solve(a == L, q)
            
                #saving solution 
                q_numr[n][:]=q.vector()[:]
                
                #Save to file 
                vtkfile_q << (q, t)
                
                #Update previous solution
                q_n.assign(q)
            
            ######## end of time stepping
        
            ## returning values from function
            return (q, q_numr)
        
        
        
        
    
        ################################################################################
        ## Function for computing q_star
        
        def q_star_computation(X, w):
        
            
            #initialize q_numr and q_star_numr 
            q_star_numr=np.zeros((n_t, n_mesh)) 
            
            ##initial q_star at t=0
            q_star_n= Function(V) 
               
            
            ##initial condition for q_star
            q_star_0= Function(V)                                      
            q_star_0_numr = (a_chain*sqrt(N_chain))*(1/(sqrt(2*pi)*c_dirac))*np.exp( (-1/(2*c_dirac**2))* (   (dof_x- X[0])**2 + (dof_y- X[1])**2  )  ) 
            q_star_0.vector()[:] = q_star_0_numr 
            
            ##initialize q_star value at t=0
            q_star_n.assign(q_star_0)
            
            #save q_star_0 in file
            vtkfile_q_star << (q_star_0, 0)
            
            #save q_star_0 in matrix
            q_star_numr[0][:]=q_star_0.vector()[:]
            
            
            
            ######## time stepping for q_star
            for n in range(1,n_t): 
                
                # print(n)
                
                t=dt*n #contor parameter value s
               
                ######################################### computing q_star
                    
                
                #defining q* and v*
                q_star = TrialFunction(V)
                v_star = TestFunction(V)
                
                #a_star and L_star in fem weak form for fenics (with Crank-Nicolson time stepping)
                a_star= G_chain*(dt/2)*dot(grad(q_star), grad(v_star))*dx + q_star*v_star*dx  + dt*w*q_star*v_star*dx 
                L_star= ( q_star_n*v_star*dx - G_chain*(dt/2)*dot(grad(q_star_n), grad(v_star))*dx  - dt*w*q_star_n*v_star*dx)
                        
                #solve variational problem
                q_star=Function(V)
                solve(a_star == L_star, q_star)
            
                #saving solution 
                q_star_numr[n][:]=q_star.vector()[:]
            
                #Save to file 
                vtkfile_q_star << (q_star, t)
                    
                #Update previous solution
                q_star_n.assign(q_star)
            
            #### end of time stepping for q_star
        
            ## returning values from function
            return (q_star, q_star_numr)
        
        
        
            
        
        
        
        ################################################################################
        ## Function for single chain computation
        
        def single_chain_computation(q_numr, q_star_numr):
               
            ##computing Q (Complete Partition Function for single chain)
            
            Q=np.zeros(n_t) # Complete Partition Function Q at each position along the chain
            
            for i in range(n_t):
                      
                q_temp1=Function(V)        
                
                q_temp1_numr=q_numr[i][:]*q_star_numr[n_t-1-i][:]
                q_temp1.vector()[:]= q_temp1_numr
                
                Q[i]=assemble(q_temp1*dx)/V_domain #Q is normalized with dividing by volume of the domain
                
            Q_chain= Q[round(n_t/2)] #Q at s=0.5  
            
            ############### computing phi (Avg segment density)
            
            phi_chain=Function(V)  # phi function
            phi_chain_numr=np.zeros(n_mesh) #phi over mesh nodes(vector of nx.ny.nz size)
            
            for i in range(n_t):
                        
                phi_chain_temp1_numr= q_numr[i][:]*q_star_numr[n_t-1-i][:]            
                phi_chain_numr= phi_chain_numr + phi_chain_temp1_numr
            
            phi_chain_numr= (1/(V_domain*Q_chain))* phi_chain_numr
                
            #converting phi into fem function                        
            phi_chain.vector()[:]= phi_chain_numr        
            
            #saving phi_chain to file 
            vtkfile_phi_chain << (phi_chain)
        
            ## returning values from function
            return (Q, phi_chain, phi_chain_numr)
        
        ################################################################################
        
        
        
        
        print('stage 1 done')
        
        
        
        
        #### computation for initial guess of w
        
        
        ###############################################################################
        ## Chain computation for Block 13
        ###############################################################################
        
        ###############################################################################
        ## Chain 1 computation
       
        #random normalized external field as initial guess
        w=Function(V)
        # w_numr=np.zeros(n_mesh)
        w_numr_temp1= np.random.normal(0, 1.0, n_mesh) 
        w_numr=w_numr_temp1/ np.max(np.abs(w_numr_temp1))
        w.vector()[:]= w_numr
        
        ## Point 1 computation
        x1= -(1)* lambda_1_stretch[i1]    * (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2)   #x component of X1 vector
        y1= -(1)* (lambda_2_stretch[i2])* (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2)  #y component of X1 vector
        X1=np.array([x1,y1]) #position vector of 1st chain start end
        # generating vtk files to store q_point_1
        vtkfile_q = File('q_point_1.pvd')
        q_point_1=Function(V) 
        q_point_1, q_point_1_numr = q_computation(X1, w)
               
        ## Point 5 computation
        x5= 0 #x component of X5 vector
        y5= 0 #y component of X5 vector
        X5=np.array([x5,y5]) #position vector of 1st chain start end
        # generating vtk files to store q_point_5
        vtkfile_q_star = File('q_point_5.pvd')
        q_point_5=Function(V) 
        q_point_5, q_point_5_numr = q_star_computation(X5, w)
        
                 
        # generating vtk file
        vtkfile_phi_chain = File('phi_chain_1.pvd')
        Q1, phi_chain_1, phi_chain_1_numr= single_chain_computation( q_point_1_numr, q_point_5_numr)
        
        ###############################################################################
        ## Chain 2 computation
       
        #random normalized external field as initial guess
        w=Function(V)
        # w_numr=np.zeros(n_mesh)
        w_numr_temp1= np.random.normal(0, 1.0, n_mesh) 
        w_numr=w_numr_temp1/ np.max(np.abs(w_numr_temp1))
        w.vector()[:]= w_numr
        
        ## Point 2 computation
        x2=  (1)* (lambda_1_stretch[i1])   * (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2)   #x component of X2 vector
        y2= -(1)* (lambda_2_stretch[i2]) * (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2)   #y component of X2 vector
        X2=np.array([x2,y2]) #position vector of 1st chain start end
        # generating vtk files to store q_point_1
        vtkfile_q = File('q_point_2.pvd')
        q_point_2=Function(V) 
        q_point_2, q_point_2_numr = q_computation(X2, w)
               
        ## Point 5 computation
        x5= 0 #x component of X5 vector
        y5= 0 #y component of X5 vector
        X5=np.array([x5,y5]) #position vector of 1st chain start end
        # generating vtk files to store q_point_5
        vtkfile_q_star = File('q_point_5.pvd')
        q_point_5=Function(V) 
        q_point_5, q_point_5_numr = q_star_computation(X5, w)
        
                 
        # generating vtk file
        vtkfile_phi_chain = File('phi_chain_2.pvd')
        Q2, phi_chain_2, phi_chain_2_numr= single_chain_computation( q_point_2_numr, q_point_5_numr)
        
        ###############################################################################
        ## Chain 3 computation
       
        #random normalized external field as initial guess
        w=Function(V)
        # w_numr=np.zeros(n_mesh)
        w_numr_temp1= np.random.normal(0, 1.0, n_mesh) 
        w_numr=w_numr_temp1/ np.max(np.abs(w_numr_temp1))
        w.vector()[:]= w_numr
        
        ## Point 3 computation
        x3=  (1)* (lambda_1_stretch[i1])    * (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2)  #x component of X1 vector
        y3=  (1)* (lambda_2_stretch[i2])  * (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) #y component of X1 vector
        X3=np.array([x3,y3]) #position vector of 1st chain start end
        # generating vtk files to store q_point_1
        vtkfile_q = File('q_point_3.pvd')
        q_point_3=Function(V) 
        q_point_3, q_point_3_numr = q_computation(X3, w)
               
        ## Point 5 computation
        x5= 0 #x component of X5 vector
        y5= 0 #y component of X5 vector
        X5=np.array([x5,y5]) #position vector of 1st chain start end
        # generating vtk files to store q_point_5
        vtkfile_q_star = File('q_point_5.pvd')
        q_point_5=Function(V) 
        q_point_5, q_point_5_numr = q_star_computation(X5, w)
        
                 
        # generating vtk file
        vtkfile_phi_chain = File('phi_chain_3.pvd')
        Q3, phi_chain_3, phi_chain_3_numr= single_chain_computation( q_point_3_numr, q_point_5_numr)
        
        ###############################################################################
        ## Chain 4 computation
       
        #random normalized external field as initial guess
        w=Function(V)
        # w_numr=np.zeros(n_mesh)
        w_numr_temp1= np.random.normal(0, 1.0, n_mesh) 
        w_numr=w_numr_temp1/ np.max(np.abs(w_numr_temp1))
        w.vector()[:]= w_numr
        
        ## Point 4 computation
        x4= -(1)* (lambda_1_stretch[i1])   * (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2)  #x component of X1 vector
        y4=  (1)* (lambda_2_stretch[i2]) * (1/(2*sqrt(2)))*a_chain*(N_chain)**(1/2) #y component of X1 vector
        X4=np.array([x4,y4]) #position vector of 1st chain start end
        # generating vtk files to store q_point_1
        vtkfile_q = File('q_point_4.pvd')
        q_point_4=Function(V) 
        q_point_4, q_point_4_numr = q_computation(X4, w)
               
        ## Point 5 computation
        x5= 0 #x component of X5 vector
        y5= 0 #y component of X5 vector
        X5=np.array([x5,y5]) #position vector of 1st chain start end
        # generating vtk files to store q_point_5
        vtkfile_q_star = File('q_point_5.pvd')
        q_point_5=Function(V) 
        q_point_5, q_point_5_numr = q_star_computation(X5, w)
        
                 
        # generating vtk file
        vtkfile_phi_chain = File('phi_chain_4.pvd')
        Q4, phi_chain_4, phi_chain_4_numr= single_chain_computation( q_point_4_numr, q_point_5_numr)
        
    
        
        
        ###############################################################################
        #getting Q at each chain mid point        
      
        Q1_mid= Q1[round(n_t/2)] 
        Q2_mid= Q2[round(n_t/2)] 
        Q3_mid= Q3[round(n_t/2)] 
        Q4_mid= Q4[round(n_t/2)]   
      
        
        
        #computing phi_MF and converting it into fem function
        
        phi_MF=Function(V)  # phi function
        phi_MF_numr=np.zeros(n_mesh) #phi over mesh nodes(vector of nx.ny.nz size)
        
        phi_MF_numr= phi_chain_1_numr + phi_chain_2_numr + phi_chain_3_numr + phi_chain_4_numr  
        phi_MF.vector()[:]= phi_MF_numr        
      
      
        # saving phi_MF in file
        vtkfile_phi_MF = File('phi_MF.pvd')
        vtkfile_phi_MF << (phi_MF)
        
        #total free energy, for external field w
        H = k_B*Temp*(1/(2*u0))*assemble(w*w*dx)  - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  )
        
        H_entropic= - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  )
        H_interaction= k_B*Temp*(1/(2*u0))*assemble(w*w*dx)
    
        
        
    
        
    
        print('stage 2 done')
    
    
    
    
        
        ############################################################### 
        ## Shifting functions 
        ############################################################### 
        
        phi_MF_B13= interpolate (phi_MF, V_B13)
        File('phi_MF_B13.pvd') << (phi_MF_B13)
    
    
        phi_MF_B1= Function(V_B1)
        phi_MF_B1.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B1.pvd') << (phi_MF_B1)
    
        phi_MF_B2= Function(V_B2)
        phi_MF_B2.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B2.pvd') << (phi_MF_B2)
    
        phi_MF_B3= Function(V_B3)
        phi_MF_B3.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B3.pvd') << (phi_MF_B3)
        
        phi_MF_B4= Function(V_B4)
        phi_MF_B4.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B4.pvd') << (phi_MF_B4)
        
        phi_MF_B5= Function(V_B5)
        phi_MF_B5.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B5.pvd') << (phi_MF_B5)
        
        phi_MF_B6= Function(V_B6)
        phi_MF_B6.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B6.pvd') << (phi_MF_B6)
        
        phi_MF_B7= Function(V_B7)
        phi_MF_B7.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B7.pvd') << (phi_MF_B7)
        
        phi_MF_B8= Function(V_B8)
        phi_MF_B8.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B8.pvd') << (phi_MF_B8)
        
        phi_MF_B9= Function(V_B9)
        phi_MF_B9.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B9.pvd') << (phi_MF_B9)
        
        phi_MF_B10= Function(V_B10)
        phi_MF_B10.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B10.pvd') << (phi_MF_B10)
        
        phi_MF_B11= Function(V_B11)
        phi_MF_B11.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B11.pvd') << (phi_MF_B11)
        
        phi_MF_B12= Function(V_B12)
        phi_MF_B12.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B12.pvd') << (phi_MF_B12)
        
        
        phi_MF_B14= Function(V_B14)
        phi_MF_B14.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B14.pvd') << (phi_MF_B14)
        
        phi_MF_B15= Function(V_B15)
        phi_MF_B15.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B15.pvd') << (phi_MF_B15)
        
        phi_MF_B16= Function(V_B16)
        phi_MF_B16.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B16.pvd') << (phi_MF_B16)
        
        phi_MF_B17= Function(V_B17)
        phi_MF_B17.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B17.pvd') << (phi_MF_B17)
        
        phi_MF_B18= Function(V_B18)
        phi_MF_B18.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B18.pvd') << (phi_MF_B18)
        
        phi_MF_B19= Function(V_B19)
        phi_MF_B19.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B19.pvd') << (phi_MF_B19)
        
        phi_MF_B20= Function(V_B20)
        phi_MF_B20.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B20.pvd') << (phi_MF_B20)
        
        phi_MF_B21= Function(V_B21)
        phi_MF_B21.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B21.pvd') << (phi_MF_B21)
        
        phi_MF_B22= Function(V_B22)
        phi_MF_B22.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B22.pvd') << (phi_MF_B22)
        
        phi_MF_B23= Function(V_B23)
        phi_MF_B23.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B23.pvd') << (phi_MF_B23)
        
        phi_MF_B24= Function(V_B24)
        phi_MF_B24.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B24.pvd') << (phi_MF_B24)
        
        phi_MF_B25= Function(V_B25)
        phi_MF_B25.vector()[:]= phi_MF_B13.vector()[:]
        File('phi_MF_B25.pvd') << (phi_MF_B25)
    
    
        
    
    
        print('stage 3 done')
    
    
    
    
    
        ############################################################### 
        ## Getting solutions on submesh 
        
        phi_MF_B13_submesh= interpolate (phi_MF_B13, V_B13_submesh)
        File('phi_MF_B13_submesh.pvd') << (phi_MF_B13_submesh)
    
    
        phi_MF_B1_submesh= interpolate (phi_MF_B1, V_B1_submesh)
        File('phi_MF_B1_submesh.pvd') << (phi_MF_B1_submesh)
    
    
        phi_MF_B2_submesh= interpolate (phi_MF_B2, V_B2_submesh)
        File('phi_MF_B2_submesh.pvd') << (phi_MF_B2_submesh)
        
        
        phi_MF_B3_submesh= interpolate (phi_MF_B3, V_B3_submesh)
        File('phi_MF_B3_submesh.pvd') << (phi_MF_B3_submesh)
        
        
        phi_MF_B4_submesh= interpolate (phi_MF_B4, V_B4_submesh)
        File('phi_MF_B4_submesh.pvd') << (phi_MF_B4_submesh)
        
        
        phi_MF_B5_submesh= interpolate (phi_MF_B5, V_B5_submesh)
        File('phi_MF_B5_submesh.pvd') << (phi_MF_B5_submesh)
        
        
        phi_MF_B6_submesh= interpolate (phi_MF_B6, V_B6_submesh)
        File('phi_MF_B6_submesh.pvd') << (phi_MF_B6_submesh)
        
        
        phi_MF_B7_submesh= interpolate (phi_MF_B7, V_B7_submesh)
        File('phi_MF_B7_submesh.pvd') << (phi_MF_B7_submesh)
        
        
        phi_MF_B8_submesh= interpolate (phi_MF_B8, V_B8_submesh)
        File('phi_MF_B8_submesh.pvd') << (phi_MF_B8_submesh)
        
        
        phi_MF_B9_submesh= interpolate (phi_MF_B9, V_B9_submesh)
        File('phi_MF_B9_submesh.pvd') << (phi_MF_B9_submesh)
        
        
        phi_MF_B10_submesh= interpolate (phi_MF_B10, V_B10_submesh)
        File('phi_MF_B10_submesh.pvd') << (phi_MF_B10_submesh)
        
        
        phi_MF_B11_submesh= interpolate (phi_MF_B11, V_B11_submesh)
        File('phi_MF_B11_submesh.pvd') << (phi_MF_B11_submesh)
        
        
        phi_MF_B12_submesh= interpolate (phi_MF_B12, V_B12_submesh)
        File('phi_MF_B12_submesh.pvd') << (phi_MF_B12_submesh)
        
        
        
        phi_MF_B14_submesh= interpolate (phi_MF_B14, V_B14_submesh)
        File('phi_MF_B14_submesh.pvd') << (phi_MF_B14_submesh)
        
        
        phi_MF_B15_submesh= interpolate (phi_MF_B15, V_B15_submesh)
        File('phi_MF_B15_submesh.pvd') << (phi_MF_B15_submesh)
        
        
        phi_MF_B16_submesh= interpolate (phi_MF_B16, V_B16_submesh)
        File('phi_MF_B16_submesh.pvd') << (phi_MF_B16_submesh)
        
        
        phi_MF_B17_submesh= interpolate (phi_MF_B17, V_B17_submesh)
        File('phi_MF_B17_submesh.pvd') << (phi_MF_B17_submesh)
        
        
        phi_MF_B18_submesh= interpolate (phi_MF_B18, V_B18_submesh)
        File('phi_MF_B18_submesh.pvd') << (phi_MF_B18_submesh)
        
        
        phi_MF_B19_submesh= interpolate (phi_MF_B19, V_B19_submesh)
        File('phi_MF_B19_submesh.pvd') << (phi_MF_B19_submesh)
        
        
        phi_MF_B20_submesh= interpolate (phi_MF_B20, V_B20_submesh)
        File('phi_MF_B20_submesh.pvd') << (phi_MF_B20_submesh)
        
        
        phi_MF_B21_submesh= interpolate (phi_MF_B21, V_B21_submesh)
        File('phi_MF_B21_submesh.pvd') << (phi_MF_B21_submesh)
        
        
        phi_MF_B22_submesh= interpolate (phi_MF_B22, V_B22_submesh)
        File('phi_MF_B22_submesh.pvd') << (phi_MF_B22_submesh)
        
        
        phi_MF_B23_submesh= interpolate (phi_MF_B23, V_B23_submesh)
        File('phi_MF_B23_submesh.pvd') << (phi_MF_B23_submesh)
        
        
        phi_MF_B24_submesh= interpolate (phi_MF_B24, V_B24_submesh)
        File('phi_MF_B24_submesh.pvd') << (phi_MF_B24_submesh)
        
        
        phi_MF_B25_submesh= interpolate (phi_MF_B25, V_B25_submesh)
        File('phi_MF_B25_submesh.pvd') << (phi_MF_B25_submesh)
    
    
    
    
    
        print('stage 4 done')
    
    
    
    
    
        ####################################### converting submesh functions to big-box-mesh function    
        
        
        ##################### getting phi_MF_B13_submesh on V_big
    
        dof_v_V_B13_submesh= np.array(dof_to_vertex_map(V_B13_submesh), dtype=int) 
        v_V_B13_submesh = mesh_B13_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B13_submesh_on_V_big = Function(V_big)
        array = phi_MF_B13_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B13_submesh[dof_v_V_B13_submesh]]] = phi_MF_B13_submesh.vector()[:]
        phi_MF_B13_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B13_submesh_on_V_big.pvd') << (phi_MF_B13_submesh_on_V_big)
    
        
        ##################### getting phi_MF_B1_submesh on V_big
    
        dof_v_V_B1_submesh= np.array(dof_to_vertex_map(V_B1_submesh), dtype=int) 
        v_V_B1_submesh = mesh_B1_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B1_submesh_on_V_big = Function(V_big)
        array = phi_MF_B1_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B1_submesh[dof_v_V_B1_submesh]]] = phi_MF_B1_submesh.vector()[:]
        phi_MF_B1_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B1_submesh_on_V_big.pvd') << (phi_MF_B1_submesh_on_V_big)
        
        ##################### getting phi_MF_B2_submesh on V_big
    
        dof_v_V_B2_submesh = np.array(dof_to_vertex_map(V_B2_submesh), dtype=int) 
        v_V_B2_submesh = mesh_B2_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B2_submesh_on_V_big = Function(V_big)
        array = phi_MF_B2_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B2_submesh[dof_v_V_B2_submesh]]] = phi_MF_B2_submesh.vector()[:]
        phi_MF_B2_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B2_submesh_on_V_big.pvd') << (phi_MF_B2_submesh_on_V_big)
        
        ##################### getting phi_MF_B3_submesh on V_big
    
        dof_v_V_B3_submesh= np.array(dof_to_vertex_map(V_B3_submesh), dtype=int) 
        v_V_B3_submesh = mesh_B3_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B3_submesh_on_V_big = Function(V_big)
        array = phi_MF_B3_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B3_submesh[dof_v_V_B3_submesh]]] = phi_MF_B3_submesh.vector()[:]
        phi_MF_B3_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B3_submesh_on_V_big.pvd') << (phi_MF_B3_submesh_on_V_big)
        
        ##################### getting phi_MF_B4_submesh on V_big
    
        dof_v_V_B4_submesh= np.array(dof_to_vertex_map(V_B4_submesh), dtype=int) 
        v_V_B4_submesh = mesh_B4_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B4_submesh_on_V_big = Function(V_big)
        array = phi_MF_B4_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B4_submesh[dof_v_V_B4_submesh]]] = phi_MF_B4_submesh.vector()[:]
        phi_MF_B4_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B4_submesh_on_V_big.pvd') << (phi_MF_B4_submesh_on_V_big)
        
        ##################### getting phi_MF_B5_submesh on V_big
    
        dof_v_V_B5_submesh= np.array(dof_to_vertex_map(V_B5_submesh), dtype=int) 
        v_V_B5_submesh = mesh_B5_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B5_submesh_on_V_big = Function(V_big)
        array = phi_MF_B5_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B5_submesh[dof_v_V_B5_submesh]]] = phi_MF_B5_submesh.vector()[:]
        phi_MF_B5_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B5_submesh_on_V_big.pvd') << (phi_MF_B5_submesh_on_V_big)
        
        ##################### getting phi_MF_B6_submesh on V_big
    
        dof_v_V_B6_submesh= np.array(dof_to_vertex_map(V_B6_submesh), dtype=int) 
        v_V_B6_submesh = mesh_B6_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B6_submesh_on_V_big = Function(V_big)
        array = phi_MF_B6_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B6_submesh[dof_v_V_B6_submesh]]] = phi_MF_B6_submesh.vector()[:]
        phi_MF_B6_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B6_submesh_on_V_big.pvd') << (phi_MF_B6_submesh_on_V_big)
        
        ##################### getting phi_MF_B7_submesh on V_big
    
        dof_v_V_B7_submesh= np.array(dof_to_vertex_map(V_B7_submesh), dtype=int) 
        v_V_B7_submesh = mesh_B7_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B7_submesh_on_V_big = Function(V_big)
        array = phi_MF_B7_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B7_submesh[dof_v_V_B7_submesh]]] = phi_MF_B7_submesh.vector()[:]
        phi_MF_B7_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B7_submesh_on_V_big.pvd') << (phi_MF_B7_submesh_on_V_big)
        
        ##################### getting phi_MF_B8_submesh on V_big
    
        dof_v_V_B8_submesh= np.array(dof_to_vertex_map(V_B8_submesh), dtype=int) 
        v_V_B8_submesh = mesh_B8_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B8_submesh_on_V_big = Function(V_big)
        array = phi_MF_B8_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B8_submesh[dof_v_V_B8_submesh]]] = phi_MF_B8_submesh.vector()[:]
        phi_MF_B8_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B8_submesh_on_V_big.pvd') << (phi_MF_B8_submesh_on_V_big)
        
        ##################### getting phi_MF_B9_submesh on V_big
    
        dof_v_V_B9_submesh= np.array(dof_to_vertex_map(V_B9_submesh), dtype=int) 
        v_V_B9_submesh = mesh_B9_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B9_submesh_on_V_big = Function(V_big)
        array = phi_MF_B9_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B9_submesh[dof_v_V_B9_submesh]]] = phi_MF_B9_submesh.vector()[:]
        phi_MF_B9_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B9_submesh_on_V_big.pvd') << (phi_MF_B9_submesh_on_V_big)
        
        ##################### getting phi_MF_B10_submesh on V_big
    
        dof_v_V_B10_submesh= np.array(dof_to_vertex_map(V_B10_submesh), dtype=int) 
        v_V_B10_submesh = mesh_B10_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B10_submesh_on_V_big = Function(V_big)
        array = phi_MF_B10_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B10_submesh[dof_v_V_B10_submesh]]] = phi_MF_B10_submesh.vector()[:]
        phi_MF_B10_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B10_submesh_on_V_big.pvd') << (phi_MF_B10_submesh_on_V_big)
        
        ##################### getting phi_MF_B11_submesh on V_big
    
        dof_v_V_B11_submesh= np.array(dof_to_vertex_map(V_B11_submesh), dtype=int) 
        v_V_B11_submesh = mesh_B11_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B11_submesh_on_V_big = Function(V_big)
        array = phi_MF_B11_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B11_submesh[dof_v_V_B11_submesh]]] = phi_MF_B11_submesh.vector()[:]
        phi_MF_B11_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B11_submesh_on_V_big.pvd') << (phi_MF_B11_submesh_on_V_big)
        
        ##################### getting phi_MF_B12_submesh on V_big
    
        dof_v_V_B12_submesh= np.array(dof_to_vertex_map(V_B12_submesh), dtype=int) 
        v_V_B12_submesh = mesh_B12_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B12_submesh_on_V_big = Function(V_big)
        array = phi_MF_B12_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B12_submesh[dof_v_V_B12_submesh]]] = phi_MF_B12_submesh.vector()[:]
        phi_MF_B12_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B12_submesh_on_V_big.pvd') << (phi_MF_B12_submesh_on_V_big)
        
        ##################### getting phi_MF_B13_submesh on V_big
    
        dof_v_V_B13_submesh= np.array(dof_to_vertex_map(V_B13_submesh), dtype=int) 
        v_V_B13_submesh = mesh_B13_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B13_submesh_on_V_big = Function(V_big)
        array = phi_MF_B13_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B13_submesh[dof_v_V_B13_submesh]]] = phi_MF_B13_submesh.vector()[:]
        phi_MF_B13_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B13_submesh_on_V_big.pvd') << (phi_MF_B13_submesh_on_V_big)
        
        ##################### getting phi_MF_B14_submesh on V_big
    
        dof_v_V_B14_submesh= np.array(dof_to_vertex_map(V_B14_submesh), dtype=int) 
        v_V_B14_submesh = mesh_B14_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B14_submesh_on_V_big = Function(V_big)
        array = phi_MF_B14_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B14_submesh[dof_v_V_B14_submesh]]] = phi_MF_B14_submesh.vector()[:]
        phi_MF_B14_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B14_submesh_on_V_big.pvd') << (phi_MF_B14_submesh_on_V_big)
        
        ##################### getting phi_MF_B15_submesh on V_big
    
        dof_v_V_B15_submesh= np.array(dof_to_vertex_map(V_B15_submesh), dtype=int) 
        v_V_B15_submesh = mesh_B15_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B15_submesh_on_V_big = Function(V_big)
        array = phi_MF_B15_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B15_submesh[dof_v_V_B15_submesh]]] = phi_MF_B15_submesh.vector()[:]
        phi_MF_B15_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B15_submesh_on_V_big.pvd') << (phi_MF_B15_submesh_on_V_big)
        
        ##################### getting phi_MF_B16_submesh on V_big
    
        dof_v_V_B16_submesh= np.array(dof_to_vertex_map(V_B16_submesh), dtype=int) 
        v_V_B16_submesh = mesh_B16_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B16_submesh_on_V_big = Function(V_big)
        array = phi_MF_B16_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B16_submesh[dof_v_V_B16_submesh]]] = phi_MF_B16_submesh.vector()[:]
        phi_MF_B16_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B16_submesh_on_V_big.pvd') << (phi_MF_B16_submesh_on_V_big)
        
        ##################### getting phi_MF_B17_submesh on V_big
    
        dof_v_V_B17_submesh= np.array(dof_to_vertex_map(V_B17_submesh), dtype=int) 
        v_V_B17_submesh = mesh_B17_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B17_submesh_on_V_big = Function(V_big)
        array = phi_MF_B17_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B17_submesh[dof_v_V_B17_submesh]]] = phi_MF_B17_submesh.vector()[:]
        phi_MF_B17_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B17_submesh_on_V_big.pvd') << (phi_MF_B17_submesh_on_V_big)
        
        ##################### getting phi_MF_B18_submesh on V_big
    
        dof_v_V_B18_submesh= np.array(dof_to_vertex_map(V_B18_submesh), dtype=int) 
        v_V_B18_submesh = mesh_B18_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B18_submesh_on_V_big = Function(V_big)
        array = phi_MF_B18_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B18_submesh[dof_v_V_B18_submesh]]] = phi_MF_B18_submesh.vector()[:]
        phi_MF_B18_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B18_submesh_on_V_big.pvd') << (phi_MF_B18_submesh_on_V_big)
        
        ##################### getting phi_MF_B19_submesh on V_big
    
        dof_v_V_B19_submesh= np.array(dof_to_vertex_map(V_B19_submesh), dtype=int) 
        v_V_B19_submesh = mesh_B19_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B19_submesh_on_V_big = Function(V_big)
        array = phi_MF_B19_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B19_submesh[dof_v_V_B19_submesh]]] = phi_MF_B19_submesh.vector()[:]
        phi_MF_B19_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B19_submesh_on_V_big.pvd') << (phi_MF_B19_submesh_on_V_big)
        
        ##################### getting phi_MF_B20_submesh on V_big
    
        dof_v_V_B20_submesh= np.array(dof_to_vertex_map(V_B20_submesh), dtype=int) 
        v_V_B20_submesh = mesh_B20_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B20_submesh_on_V_big = Function(V_big)
        array = phi_MF_B20_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B20_submesh[dof_v_V_B20_submesh]]] = phi_MF_B20_submesh.vector()[:]
        phi_MF_B20_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B20_submesh_on_V_big.pvd') << (phi_MF_B20_submesh_on_V_big)
        
        ##################### getting phi_MF_B21_submesh on V_big
    
        dof_v_V_B21_submesh= np.array(dof_to_vertex_map(V_B21_submesh), dtype=int) 
        v_V_B21_submesh = mesh_B21_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B21_submesh_on_V_big = Function(V_big)
        array = phi_MF_B21_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B21_submesh[dof_v_V_B21_submesh]]] = phi_MF_B21_submesh.vector()[:]
        phi_MF_B21_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B21_submesh_on_V_big.pvd') << (phi_MF_B21_submesh_on_V_big)
        
        ##################### getting phi_MF_B22_submesh on V_big
    
        dof_v_V_B22_submesh= np.array(dof_to_vertex_map(V_B22_submesh), dtype=int) 
        v_V_B22_submesh = mesh_B22_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B22_submesh_on_V_big = Function(V_big)
        array = phi_MF_B22_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B22_submesh[dof_v_V_B22_submesh]]] = phi_MF_B22_submesh.vector()[:]
        phi_MF_B22_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B22_submesh_on_V_big.pvd') << (phi_MF_B22_submesh_on_V_big)
        
        ##################### getting phi_MF_B23_submesh on V_big
    
        dof_v_V_B23_submesh= np.array(dof_to_vertex_map(V_B23_submesh), dtype=int) 
        v_V_B23_submesh = mesh_B23_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B23_submesh_on_V_big = Function(V_big)
        array = phi_MF_B23_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B23_submesh[dof_v_V_B23_submesh]]] = phi_MF_B23_submesh.vector()[:]
        phi_MF_B23_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B23_submesh_on_V_big.pvd') << (phi_MF_B23_submesh_on_V_big)
        
        ##################### getting phi_MF_B24_submesh on V_big
    
        dof_v_V_B24_submesh= np.array(dof_to_vertex_map(V_B24_submesh), dtype=int) 
        v_V_B24_submesh = mesh_B24_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B24_submesh_on_V_big = Function(V_big)
        array = phi_MF_B24_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B24_submesh[dof_v_V_B24_submesh]]] = phi_MF_B24_submesh.vector()[:]
        phi_MF_B24_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B24_submesh_on_V_big.pvd') << (phi_MF_B24_submesh_on_V_big)
        
        ##################### getting phi_MF_B25_submesh on V_big
    
        dof_v_V_B25_submesh= np.array(dof_to_vertex_map(V_B25_submesh), dtype=int) 
        v_V_B25_submesh = mesh_B25_submesh.data().array('parent_vertex_indices', 0)
        
        v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
        
        phi_MF_B25_submesh_on_V_big = Function(V_big)
        array = phi_MF_B25_submesh_on_V_big.vector()[:]
        # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
        array[v_dof_V_big[v_V_B25_submesh[dof_v_V_B25_submesh]]] = phi_MF_B25_submesh.vector()[:]
        phi_MF_B25_submesh_on_V_big.vector()[:] = array
        
        File('phi_MF_B25_submesh_on_V_big.pvd') << (phi_MF_B25_submesh_on_V_big)
    
    
    
    
    
        print('stage 5 done')
    
    
    
        
        
        ############################################################### 
        ## getting w for next time step
        ############################################################### 
        
        phi_MF_on_V_big= Function(V_big)
        phi_MF_on_V_big.vector()[:] = phi_MF_B1_submesh_on_V_big.vector()[:] +  phi_MF_B2_submesh_on_V_big.vector()[:] +  phi_MF_B3_submesh_on_V_big.vector()[:] +  phi_MF_B4_submesh_on_V_big.vector()[:] +  phi_MF_B5_submesh_on_V_big.vector()[:] + phi_MF_B6_submesh_on_V_big.vector()[:] +  phi_MF_B7_submesh_on_V_big.vector()[:] + phi_MF_B8_submesh_on_V_big.vector()[:] +  phi_MF_B9_submesh_on_V_big.vector()[:]  + phi_MF_B10_submesh_on_V_big.vector()[:] + phi_MF_B11_submesh_on_V_big.vector()[:] + phi_MF_B12_submesh_on_V_big.vector()[:] + phi_MF_B13_submesh_on_V_big.vector()[:] + phi_MF_B14_submesh_on_V_big.vector()[:] + phi_MF_B15_submesh_on_V_big.vector()[:] + phi_MF_B16_submesh_on_V_big.vector()[:] + phi_MF_B17_submesh_on_V_big.vector()[:] + phi_MF_B18_submesh_on_V_big.vector()[:] + phi_MF_B19_submesh_on_V_big.vector()[:] + phi_MF_B20_submesh_on_V_big.vector()[:] + phi_MF_B21_submesh_on_V_big.vector()[:] + phi_MF_B22_submesh_on_V_big.vector()[:] + phi_MF_B23_submesh_on_V_big.vector()[:] + phi_MF_B24_submesh_on_V_big.vector()[:] + phi_MF_B25_submesh_on_V_big.vector()[:]  
    
        
        #saving phi_MF 
        vtkfile_phi_MF = File('phi_MF_on_V_big.pvd')
        vtkfile_phi_MF << (phi_MF_on_V_big)
        
        ########## getting phi_MF on V
        phi_MF_on_V_big_to_V= interpolate (phi_MF_on_V_big, V)
        File('phi_MF_on_V_big_to_V.pvd') << (phi_MF_on_V_big_to_V)
        
        
        #defining w field for next step
        w=Function(V)  # phi function
        w.vector()[:]= u0* phi_MF_on_V_big_to_V.vector()[:] 
    
    
    

    
    
        ###############################################################  Iterating for finding equilibrium mean field w
        count=0
        # while error_w_norm > error_w_norm_threshold:
        while abs(delta_H_ratio) > delta_H_ratio_threshold:
            count=count+1
            print(count)
            print(delta_H_ratio)
            
        
            if count == 61:
                print('count=61')
                break
        
            ###############################################################################
            ## Point 1 computation
            
            # generating vtk files to store q_point_1
            vtkfile_q = File('q_point_1.pvd')
            q_point_1=Function(V) 
            q_point_1, q_point_1_numr = q_computation(X1, w)
            
            ###############################################################################
            ## Point 2 computation
            
            # generating vtk files to store q_point_2
            vtkfile_q = File('q_point_2.pvd')
            q_point_2=Function(V)         
            q_point_2, q_point_2_numr = q_computation(X2, w)
        
            ###############################################################################
            ## Point 3 computation
            
            # generating vtk files to store q_point_3
            vtkfile_q = File('q_point_3.pvd')
            q_point_3=Function(V) 
            q_point_3, q_point_3_numr = q_computation(X3, w)
        
            ###############################################################################
            ## Point 4 computation
            
            # generating vtk files to store q_point_4
            vtkfile_q_star = File('q_point_4.pvd')
            q_point_4=Function(V) 
            q_point_4, q_point_4_numr = q_computation(X4, w)
            
            ###############################################################################
            ## Point 5 computation
            
            # generating vtk files to store q_point_5
            vtkfile_q_star = File('q_point_5.pvd')
            q_point_5=Function(V) 
            q_point_5, q_point_5_numr = q_star_computation(X5, w)
        
        
        
        
            ###############################################################################
            ## Chain computation 
            ###############################################################################
        
        
            ###############################################################################
            ## Chain 1 computation
                    
            # generating vtk files to store q_point_1
            vtkfile_phi_chain = File('phi_chain_1.pvd')
           
            Q1, phi_chain_1, phi_chain_1_numr= single_chain_computation( q_point_1_numr, q_point_5_numr)
            
            ###############################################################################
            ## Chain 2 computation
            
            # generating vtk files to store q_point_2
            vtkfile_phi_chain = File('phi_chain_2.pvd')
            
            Q2, phi_chain_2, phi_chain_2_numr= single_chain_computation( q_point_2_numr, q_point_5_numr)
            
            ###############################################################################
            ## Chain 3 computation
            
            # generating vtk files to store q_point_1
            vtkfile_phi_chain = File('phi_chain_3.pvd')
            
            Q3, phi_chain_3, phi_chain_3_numr= single_chain_computation( q_point_3_numr, q_point_5_numr)        
            
            ###############################################################################
            ## Chain 4 computation
            
            # generating vtk files to store q_point_4
            vtkfile_phi_chain = File('phi_chain_4.pvd')
            
            Q4, phi_chain_4, phi_chain_4_numr= single_chain_computation( q_point_4_numr, q_point_5_numr)   
            
            
            
            
            
            ###############################################################################
            ##getting Q at each chain mid point        
           
            Q1_mid= Q1[round(n_t/2)] 
            Q2_mid= Q2[round(n_t/2)] 
            Q3_mid= Q3[round(n_t/2)] 
            Q4_mid= Q4[round(n_t/2)] 
           
            
          
            #computing phi_MF and converting it into fem function
            phi_MF=Function(V)  # phi function
            phi_MF_numr=np.zeros(n_mesh) #phi over mesh nodes(vector of nx.ny.nz size)
            
            phi_MF_numr= phi_chain_1_numr + phi_chain_2_numr + phi_chain_3_numr + phi_chain_4_numr                 
            phi_MF.vector()[:]= phi_MF_numr        
            
            # saving phi_MF in file
            vtkfile_phi_MF = File('phi_MF.pvd')
            vtkfile_phi_MF << (phi_MF)
            
            #total free energy, for external field w
            H = k_B*Temp*(1/(2*u0))*assemble(w*w*dx) - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  )
            H_entropic= - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  )
            H_interaction= k_B*Temp*(1/(2*u0))*assemble(w*w*dx)
            

            
            
        
        
            print('stage 2 done')
        
        
        
        
            ############################################################### 
            ## To get w for next time step
            ###############################################################
        
            
            ############################################################### 
            ## Shifting functions 
            ############################################################### 
            
            phi_MF_B13= interpolate (phi_MF, V_B13)
            File('phi_MF_B13.pvd') << (phi_MF_B13)
        
        
            phi_MF_B1= Function(V_B1)
            phi_MF_B1.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B1.pvd') << (phi_MF_B1)
        
            phi_MF_B2= Function(V_B2)
            phi_MF_B2.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B2.pvd') << (phi_MF_B2)
        
            phi_MF_B3= Function(V_B3)
            phi_MF_B3.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B3.pvd') << (phi_MF_B3)
            
            phi_MF_B4= Function(V_B4)
            phi_MF_B4.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B4.pvd') << (phi_MF_B4)
            
            phi_MF_B5= Function(V_B5)
            phi_MF_B5.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B5.pvd') << (phi_MF_B5)
            
            phi_MF_B6= Function(V_B6)
            phi_MF_B6.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B6.pvd') << (phi_MF_B6)
            
            phi_MF_B7= Function(V_B7)
            phi_MF_B7.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B7.pvd') << (phi_MF_B7)
            
            phi_MF_B8= Function(V_B8)
            phi_MF_B8.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B8.pvd') << (phi_MF_B8)
            
            phi_MF_B9= Function(V_B9)
            phi_MF_B9.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B9.pvd') << (phi_MF_B9)
            
            phi_MF_B10= Function(V_B10)
            phi_MF_B10.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B10.pvd') << (phi_MF_B10)
            
            phi_MF_B11= Function(V_B11)
            phi_MF_B11.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B11.pvd') << (phi_MF_B11)
            
            phi_MF_B12= Function(V_B12)
            phi_MF_B12.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B12.pvd') << (phi_MF_B12)
            
            
            phi_MF_B14= Function(V_B14)
            phi_MF_B14.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B14.pvd') << (phi_MF_B14)
            
            phi_MF_B15= Function(V_B15)
            phi_MF_B15.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B15.pvd') << (phi_MF_B15)
            
            phi_MF_B16= Function(V_B16)
            phi_MF_B16.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B16.pvd') << (phi_MF_B16)
            
            phi_MF_B17= Function(V_B17)
            phi_MF_B17.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B17.pvd') << (phi_MF_B17)
            
            phi_MF_B18= Function(V_B18)
            phi_MF_B18.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B18.pvd') << (phi_MF_B18)
            
            phi_MF_B19= Function(V_B19)
            phi_MF_B19.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B19.pvd') << (phi_MF_B19)
            
            phi_MF_B20= Function(V_B20)
            phi_MF_B20.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B20.pvd') << (phi_MF_B20)
            
            phi_MF_B21= Function(V_B21)
            phi_MF_B21.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B21.pvd') << (phi_MF_B21)
            
            phi_MF_B22= Function(V_B22)
            phi_MF_B22.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B22.pvd') << (phi_MF_B22)
            
            phi_MF_B23= Function(V_B23)
            phi_MF_B23.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B23.pvd') << (phi_MF_B23)
            
            phi_MF_B24= Function(V_B24)
            phi_MF_B24.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B24.pvd') << (phi_MF_B24)
            
            phi_MF_B25= Function(V_B25)
            phi_MF_B25.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B25.pvd') << (phi_MF_B25)
        
        
        
        
        
            print('stage 3 done')
        
        
        
        
        
            ############################################################### 
            ## Getting solutions on submesh 
            
            phi_MF_B13_submesh= interpolate (phi_MF_B13, V_B13_submesh)
            File('phi_MF_B13_submesh.pvd') << (phi_MF_B13_submesh)
        
        
            phi_MF_B1_submesh= interpolate (phi_MF_B1, V_B1_submesh)
            File('phi_MF_B1_submesh.pvd') << (phi_MF_B1_submesh)
        
        
            phi_MF_B2_submesh= interpolate (phi_MF_B2, V_B2_submesh)
            File('phi_MF_B2_submesh.pvd') << (phi_MF_B2_submesh)
            
            
            phi_MF_B3_submesh= interpolate (phi_MF_B3, V_B3_submesh)
            File('phi_MF_B3_submesh.pvd') << (phi_MF_B3_submesh)
            
            
            phi_MF_B4_submesh= interpolate (phi_MF_B4, V_B4_submesh)
            File('phi_MF_B4_submesh.pvd') << (phi_MF_B4_submesh)
            
            
            phi_MF_B5_submesh= interpolate (phi_MF_B5, V_B5_submesh)
            File('phi_MF_B5_submesh.pvd') << (phi_MF_B5_submesh)
            
            
            phi_MF_B6_submesh= interpolate (phi_MF_B6, V_B6_submesh)
            File('phi_MF_B6_submesh.pvd') << (phi_MF_B6_submesh)
            
            
            phi_MF_B7_submesh= interpolate (phi_MF_B7, V_B7_submesh)
            File('phi_MF_B7_submesh.pvd') << (phi_MF_B7_submesh)
            
            
            phi_MF_B8_submesh= interpolate (phi_MF_B8, V_B8_submesh)
            File('phi_MF_B8_submesh.pvd') << (phi_MF_B8_submesh)
            
            
            phi_MF_B9_submesh= interpolate (phi_MF_B9, V_B9_submesh)
            File('phi_MF_B9_submesh.pvd') << (phi_MF_B9_submesh)
            
            
            phi_MF_B10_submesh= interpolate (phi_MF_B10, V_B10_submesh)
            File('phi_MF_B10_submesh.pvd') << (phi_MF_B10_submesh)
            
            
            phi_MF_B11_submesh= interpolate (phi_MF_B11, V_B11_submesh)
            File('phi_MF_B11_submesh.pvd') << (phi_MF_B11_submesh)
            
            
            phi_MF_B12_submesh= interpolate (phi_MF_B12, V_B12_submesh)
            File('phi_MF_B12_submesh.pvd') << (phi_MF_B12_submesh)
            
            
            
            phi_MF_B14_submesh= interpolate (phi_MF_B14, V_B14_submesh)
            File('phi_MF_B14_submesh.pvd') << (phi_MF_B14_submesh)
            
            
            phi_MF_B15_submesh= interpolate (phi_MF_B15, V_B15_submesh)
            File('phi_MF_B15_submesh.pvd') << (phi_MF_B15_submesh)
            
            
            phi_MF_B16_submesh= interpolate (phi_MF_B16, V_B16_submesh)
            File('phi_MF_B16_submesh.pvd') << (phi_MF_B16_submesh)
            
            
            phi_MF_B17_submesh= interpolate (phi_MF_B17, V_B17_submesh)
            File('phi_MF_B17_submesh.pvd') << (phi_MF_B17_submesh)
            
            
            phi_MF_B18_submesh= interpolate (phi_MF_B18, V_B18_submesh)
            File('phi_MF_B18_submesh.pvd') << (phi_MF_B18_submesh)
            
            
            phi_MF_B19_submesh= interpolate (phi_MF_B19, V_B19_submesh)
            File('phi_MF_B19_submesh.pvd') << (phi_MF_B19_submesh)
            
            
            phi_MF_B20_submesh= interpolate (phi_MF_B20, V_B20_submesh)
            File('phi_MF_B20_submesh.pvd') << (phi_MF_B20_submesh)
            
            
            phi_MF_B21_submesh= interpolate (phi_MF_B21, V_B21_submesh)
            File('phi_MF_B21_submesh.pvd') << (phi_MF_B21_submesh)
            
            
            phi_MF_B22_submesh= interpolate (phi_MF_B22, V_B22_submesh)
            File('phi_MF_B22_submesh.pvd') << (phi_MF_B22_submesh)
            
            
            phi_MF_B23_submesh= interpolate (phi_MF_B23, V_B23_submesh)
            File('phi_MF_B23_submesh.pvd') << (phi_MF_B23_submesh)
            
            
            phi_MF_B24_submesh= interpolate (phi_MF_B24, V_B24_submesh)
            File('phi_MF_B24_submesh.pvd') << (phi_MF_B24_submesh)
            
            
            phi_MF_B25_submesh= interpolate (phi_MF_B25, V_B25_submesh)
            File('phi_MF_B25_submesh.pvd') << (phi_MF_B25_submesh)
        
        
        
        
        
            print('stage 4 done')
        
        
        
        
        
            ####################################### converting submesh functions to big-box-mesh function    
            
            
            ##################### getting phi_MF_B13_submesh on V_big
        
            dof_v_V_B13_submesh= np.array(dof_to_vertex_map(V_B13_submesh), dtype=int) 
            v_V_B13_submesh = mesh_B13_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B13_submesh_on_V_big = Function(V_big)
            array = phi_MF_B13_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B13_submesh[dof_v_V_B13_submesh]]] = phi_MF_B13_submesh.vector()[:]
            phi_MF_B13_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B13_submesh_on_V_big.pvd') << (phi_MF_B13_submesh_on_V_big)
        
            
            ##################### getting phi_MF_B1_submesh on V_big
        
            dof_v_V_B1_submesh= np.array(dof_to_vertex_map(V_B1_submesh), dtype=int) 
            v_V_B1_submesh = mesh_B1_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B1_submesh_on_V_big = Function(V_big)
            array = phi_MF_B1_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B1_submesh[dof_v_V_B1_submesh]]] = phi_MF_B1_submesh.vector()[:]
            phi_MF_B1_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B1_submesh_on_V_big.pvd') << (phi_MF_B1_submesh_on_V_big)
            
            ##################### getting phi_MF_B2_submesh on V_big
        
            dof_v_V_B2_submesh = np.array(dof_to_vertex_map(V_B2_submesh), dtype=int) 
            v_V_B2_submesh = mesh_B2_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B2_submesh_on_V_big = Function(V_big)
            array = phi_MF_B2_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B2_submesh[dof_v_V_B2_submesh]]] = phi_MF_B2_submesh.vector()[:]
            phi_MF_B2_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B2_submesh_on_V_big.pvd') << (phi_MF_B2_submesh_on_V_big)
            
            ##################### getting phi_MF_B3_submesh on V_big
        
            dof_v_V_B3_submesh= np.array(dof_to_vertex_map(V_B3_submesh), dtype=int) 
            v_V_B3_submesh = mesh_B3_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B3_submesh_on_V_big = Function(V_big)
            array = phi_MF_B3_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B3_submesh[dof_v_V_B3_submesh]]] = phi_MF_B3_submesh.vector()[:]
            phi_MF_B3_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B3_submesh_on_V_big.pvd') << (phi_MF_B3_submesh_on_V_big)
            
            ##################### getting phi_MF_B4_submesh on V_big
        
            dof_v_V_B4_submesh= np.array(dof_to_vertex_map(V_B4_submesh), dtype=int) 
            v_V_B4_submesh = mesh_B4_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B4_submesh_on_V_big = Function(V_big)
            array = phi_MF_B4_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B4_submesh[dof_v_V_B4_submesh]]] = phi_MF_B4_submesh.vector()[:]
            phi_MF_B4_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B4_submesh_on_V_big.pvd') << (phi_MF_B4_submesh_on_V_big)
            
            ##################### getting phi_MF_B5_submesh on V_big
        
            dof_v_V_B5_submesh= np.array(dof_to_vertex_map(V_B5_submesh), dtype=int) 
            v_V_B5_submesh = mesh_B5_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B5_submesh_on_V_big = Function(V_big)
            array = phi_MF_B5_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B5_submesh[dof_v_V_B5_submesh]]] = phi_MF_B5_submesh.vector()[:]
            phi_MF_B5_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B5_submesh_on_V_big.pvd') << (phi_MF_B5_submesh_on_V_big)
            
            ##################### getting phi_MF_B6_submesh on V_big
        
            dof_v_V_B6_submesh= np.array(dof_to_vertex_map(V_B6_submesh), dtype=int) 
            v_V_B6_submesh = mesh_B6_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B6_submesh_on_V_big = Function(V_big)
            array = phi_MF_B6_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B6_submesh[dof_v_V_B6_submesh]]] = phi_MF_B6_submesh.vector()[:]
            phi_MF_B6_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B6_submesh_on_V_big.pvd') << (phi_MF_B6_submesh_on_V_big)
            
            ##################### getting phi_MF_B7_submesh on V_big
        
            dof_v_V_B7_submesh= np.array(dof_to_vertex_map(V_B7_submesh), dtype=int) 
            v_V_B7_submesh = mesh_B7_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B7_submesh_on_V_big = Function(V_big)
            array = phi_MF_B7_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B7_submesh[dof_v_V_B7_submesh]]] = phi_MF_B7_submesh.vector()[:]
            phi_MF_B7_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B7_submesh_on_V_big.pvd') << (phi_MF_B7_submesh_on_V_big)
            
            ##################### getting phi_MF_B8_submesh on V_big
        
            dof_v_V_B8_submesh= np.array(dof_to_vertex_map(V_B8_submesh), dtype=int) 
            v_V_B8_submesh = mesh_B8_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B8_submesh_on_V_big = Function(V_big)
            array = phi_MF_B8_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B8_submesh[dof_v_V_B8_submesh]]] = phi_MF_B8_submesh.vector()[:]
            phi_MF_B8_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B8_submesh_on_V_big.pvd') << (phi_MF_B8_submesh_on_V_big)
            
            ##################### getting phi_MF_B9_submesh on V_big
        
            dof_v_V_B9_submesh= np.array(dof_to_vertex_map(V_B9_submesh), dtype=int) 
            v_V_B9_submesh = mesh_B9_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B9_submesh_on_V_big = Function(V_big)
            array = phi_MF_B9_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B9_submesh[dof_v_V_B9_submesh]]] = phi_MF_B9_submesh.vector()[:]
            phi_MF_B9_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B9_submesh_on_V_big.pvd') << (phi_MF_B9_submesh_on_V_big)
            
            ##################### getting phi_MF_B10_submesh on V_big
        
            dof_v_V_B10_submesh= np.array(dof_to_vertex_map(V_B10_submesh), dtype=int) 
            v_V_B10_submesh = mesh_B10_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B10_submesh_on_V_big = Function(V_big)
            array = phi_MF_B10_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B10_submesh[dof_v_V_B10_submesh]]] = phi_MF_B10_submesh.vector()[:]
            phi_MF_B10_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B10_submesh_on_V_big.pvd') << (phi_MF_B10_submesh_on_V_big)
            
            ##################### getting phi_MF_B11_submesh on V_big
        
            dof_v_V_B11_submesh= np.array(dof_to_vertex_map(V_B11_submesh), dtype=int) 
            v_V_B11_submesh = mesh_B11_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B11_submesh_on_V_big = Function(V_big)
            array = phi_MF_B11_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B11_submesh[dof_v_V_B11_submesh]]] = phi_MF_B11_submesh.vector()[:]
            phi_MF_B11_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B11_submesh_on_V_big.pvd') << (phi_MF_B11_submesh_on_V_big)
            
            ##################### getting phi_MF_B12_submesh on V_big
        
            dof_v_V_B12_submesh= np.array(dof_to_vertex_map(V_B12_submesh), dtype=int) 
            v_V_B12_submesh = mesh_B12_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B12_submesh_on_V_big = Function(V_big)
            array = phi_MF_B12_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B12_submesh[dof_v_V_B12_submesh]]] = phi_MF_B12_submesh.vector()[:]
            phi_MF_B12_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B12_submesh_on_V_big.pvd') << (phi_MF_B12_submesh_on_V_big)
            
            ##################### getting phi_MF_B13_submesh on V_big
        
            dof_v_V_B13_submesh= np.array(dof_to_vertex_map(V_B13_submesh), dtype=int) 
            v_V_B13_submesh = mesh_B13_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B13_submesh_on_V_big = Function(V_big)
            array = phi_MF_B13_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B13_submesh[dof_v_V_B13_submesh]]] = phi_MF_B13_submesh.vector()[:]
            phi_MF_B13_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B13_submesh_on_V_big.pvd') << (phi_MF_B13_submesh_on_V_big)
            
            ##################### getting phi_MF_B14_submesh on V_big
        
            dof_v_V_B14_submesh= np.array(dof_to_vertex_map(V_B14_submesh), dtype=int) 
            v_V_B14_submesh = mesh_B14_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B14_submesh_on_V_big = Function(V_big)
            array = phi_MF_B14_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B14_submesh[dof_v_V_B14_submesh]]] = phi_MF_B14_submesh.vector()[:]
            phi_MF_B14_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B14_submesh_on_V_big.pvd') << (phi_MF_B14_submesh_on_V_big)
            
            ##################### getting phi_MF_B15_submesh on V_big
        
            dof_v_V_B15_submesh= np.array(dof_to_vertex_map(V_B15_submesh), dtype=int) 
            v_V_B15_submesh = mesh_B15_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B15_submesh_on_V_big = Function(V_big)
            array = phi_MF_B15_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B15_submesh[dof_v_V_B15_submesh]]] = phi_MF_B15_submesh.vector()[:]
            phi_MF_B15_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B15_submesh_on_V_big.pvd') << (phi_MF_B15_submesh_on_V_big)
            
            ##################### getting phi_MF_B16_submesh on V_big
        
            dof_v_V_B16_submesh= np.array(dof_to_vertex_map(V_B16_submesh), dtype=int) 
            v_V_B16_submesh = mesh_B16_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B16_submesh_on_V_big = Function(V_big)
            array = phi_MF_B16_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B16_submesh[dof_v_V_B16_submesh]]] = phi_MF_B16_submesh.vector()[:]
            phi_MF_B16_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B16_submesh_on_V_big.pvd') << (phi_MF_B16_submesh_on_V_big)
            
            ##################### getting phi_MF_B17_submesh on V_big
        
            dof_v_V_B17_submesh= np.array(dof_to_vertex_map(V_B17_submesh), dtype=int) 
            v_V_B17_submesh = mesh_B17_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B17_submesh_on_V_big = Function(V_big)
            array = phi_MF_B17_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B17_submesh[dof_v_V_B17_submesh]]] = phi_MF_B17_submesh.vector()[:]
            phi_MF_B17_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B17_submesh_on_V_big.pvd') << (phi_MF_B17_submesh_on_V_big)
            
            ##################### getting phi_MF_B18_submesh on V_big
        
            dof_v_V_B18_submesh= np.array(dof_to_vertex_map(V_B18_submesh), dtype=int) 
            v_V_B18_submesh = mesh_B18_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B18_submesh_on_V_big = Function(V_big)
            array = phi_MF_B18_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B18_submesh[dof_v_V_B18_submesh]]] = phi_MF_B18_submesh.vector()[:]
            phi_MF_B18_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B18_submesh_on_V_big.pvd') << (phi_MF_B18_submesh_on_V_big)
            
            ##################### getting phi_MF_B19_submesh on V_big
        
            dof_v_V_B19_submesh= np.array(dof_to_vertex_map(V_B19_submesh), dtype=int) 
            v_V_B19_submesh = mesh_B19_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B19_submesh_on_V_big = Function(V_big)
            array = phi_MF_B19_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B19_submesh[dof_v_V_B19_submesh]]] = phi_MF_B19_submesh.vector()[:]
            phi_MF_B19_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B19_submesh_on_V_big.pvd') << (phi_MF_B19_submesh_on_V_big)
            
            ##################### getting phi_MF_B20_submesh on V_big
        
            dof_v_V_B20_submesh= np.array(dof_to_vertex_map(V_B20_submesh), dtype=int) 
            v_V_B20_submesh = mesh_B20_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B20_submesh_on_V_big = Function(V_big)
            array = phi_MF_B20_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B20_submesh[dof_v_V_B20_submesh]]] = phi_MF_B20_submesh.vector()[:]
            phi_MF_B20_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B20_submesh_on_V_big.pvd') << (phi_MF_B20_submesh_on_V_big)
            
            ##################### getting phi_MF_B21_submesh on V_big
        
            dof_v_V_B21_submesh= np.array(dof_to_vertex_map(V_B21_submesh), dtype=int) 
            v_V_B21_submesh = mesh_B21_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B21_submesh_on_V_big = Function(V_big)
            array = phi_MF_B21_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B21_submesh[dof_v_V_B21_submesh]]] = phi_MF_B21_submesh.vector()[:]
            phi_MF_B21_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B21_submesh_on_V_big.pvd') << (phi_MF_B21_submesh_on_V_big)
            
            ##################### getting phi_MF_B22_submesh on V_big
        
            dof_v_V_B22_submesh= np.array(dof_to_vertex_map(V_B22_submesh), dtype=int) 
            v_V_B22_submesh = mesh_B22_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B22_submesh_on_V_big = Function(V_big)
            array = phi_MF_B22_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B22_submesh[dof_v_V_B22_submesh]]] = phi_MF_B22_submesh.vector()[:]
            phi_MF_B22_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B22_submesh_on_V_big.pvd') << (phi_MF_B22_submesh_on_V_big)
            
            ##################### getting phi_MF_B23_submesh on V_big
        
            dof_v_V_B23_submesh= np.array(dof_to_vertex_map(V_B23_submesh), dtype=int) 
            v_V_B23_submesh = mesh_B23_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B23_submesh_on_V_big = Function(V_big)
            array = phi_MF_B23_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B23_submesh[dof_v_V_B23_submesh]]] = phi_MF_B23_submesh.vector()[:]
            phi_MF_B23_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B23_submesh_on_V_big.pvd') << (phi_MF_B23_submesh_on_V_big)
            
            ##################### getting phi_MF_B24_submesh on V_big
        
            dof_v_V_B24_submesh= np.array(dof_to_vertex_map(V_B24_submesh), dtype=int) 
            v_V_B24_submesh = mesh_B24_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B24_submesh_on_V_big = Function(V_big)
            array = phi_MF_B24_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B24_submesh[dof_v_V_B24_submesh]]] = phi_MF_B24_submesh.vector()[:]
            phi_MF_B24_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B24_submesh_on_V_big.pvd') << (phi_MF_B24_submesh_on_V_big)
            
            ##################### getting phi_MF_B25_submesh on V_big
        
            dof_v_V_B25_submesh= np.array(dof_to_vertex_map(V_B25_submesh), dtype=int) 
            v_V_B25_submesh = mesh_B25_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B25_submesh_on_V_big = Function(V_big)
            array = phi_MF_B25_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B25_submesh[dof_v_V_B25_submesh]]] = phi_MF_B25_submesh.vector()[:]
            phi_MF_B25_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B25_submesh_on_V_big.pvd') << (phi_MF_B25_submesh_on_V_big)
        
        
        
        
        
            print('stage 5 done')
        
        
        
            
            
            ############################################################### 
            ## getting w for next time step
            ############################################################### 
            
            phi_MF_on_V_big= Function(V_big)
            phi_MF_on_V_big.vector()[:] = phi_MF_B1_submesh_on_V_big.vector()[:] +  phi_MF_B2_submesh_on_V_big.vector()[:] +  phi_MF_B3_submesh_on_V_big.vector()[:] +  phi_MF_B4_submesh_on_V_big.vector()[:] +  phi_MF_B5_submesh_on_V_big.vector()[:] + phi_MF_B6_submesh_on_V_big.vector()[:] +  phi_MF_B7_submesh_on_V_big.vector()[:] + phi_MF_B8_submesh_on_V_big.vector()[:] +  phi_MF_B9_submesh_on_V_big.vector()[:]  + phi_MF_B10_submesh_on_V_big.vector()[:] + phi_MF_B11_submesh_on_V_big.vector()[:] + phi_MF_B12_submesh_on_V_big.vector()[:] + phi_MF_B13_submesh_on_V_big.vector()[:] + phi_MF_B14_submesh_on_V_big.vector()[:] + phi_MF_B15_submesh_on_V_big.vector()[:] + phi_MF_B16_submesh_on_V_big.vector()[:] + phi_MF_B17_submesh_on_V_big.vector()[:] + phi_MF_B18_submesh_on_V_big.vector()[:] + phi_MF_B19_submesh_on_V_big.vector()[:] + phi_MF_B20_submesh_on_V_big.vector()[:] + phi_MF_B21_submesh_on_V_big.vector()[:] + phi_MF_B22_submesh_on_V_big.vector()[:] + phi_MF_B23_submesh_on_V_big.vector()[:] + phi_MF_B24_submesh_on_V_big.vector()[:] + phi_MF_B25_submesh_on_V_big.vector()[:]  
        
            
            #saving phi_MF 
            vtkfile_phi_MF = File('phi_MF_on_V_big.pvd')
            vtkfile_phi_MF << (phi_MF_on_V_big)
            
            ########## getting phi_MF on V
            phi_MF_on_V_big_to_V= interpolate (phi_MF_on_V_big, V)
            File('phi_MF_on_V_big_to_V.pvd') << (phi_MF_on_V_big_to_V)
            
            
            
            
            
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
            w_next.vector()[:]= u0* phi_MF_on_V_big_to_V.vector()[:] 
            
            
            
            
            
            ###############################################################################
            ## q computation for w_next
            ###############################################################################
        
            
            ###############################################################################
            ## Point 1 computation 
            
            # generating vtk files to store q_point_1
            vtkfile_q = File('q_point_1.pvd')
            q_point_1, q_point_1_numr = q_computation(X1, w_next)
            
            ###############################################################################
            ## Point 2 computation  
            
            # generating vtk files to store q_point_2
            vtkfile_q = File('q_point_2.pvd')
            q_point_2, q_point_2_numr = q_computation(X2, w_next)
        
            ###############################################################################
            ## Point 3 computation       
            
            # generating vtk files to store q_point_3
            vtkfile_q = File('q_point_3.pvd')
            q_point_3, q_point_3_numr = q_computation(X3, w_next)
        
            ###############################################################################
            ## Point 4 computation       
            
            # generating vtk files to store q_point_4
            vtkfile_q = File('q_point_4.pvd')
            q_point_4, q_point_4_numr = q_computation(X4, w_next)
        
            ###############################################################################
            ## Point 5 computation       
            
            # generating vtk files to store q_point_5
            vtkfile_q = File('q_point_5.pvd')
            q_point_5, q_point_5_numr = q_star_computation(X5, w_next)
        
        
        
        
            ###############################################################################
            ## Chain computation for w_next
            ###############################################################################
        
        
            ###############################################################################
            ## Chain 1 computation
                    
            # generating vtk files to store q_point_1
            vtkfile_phi_chain = File('phi_chain_1.pvd')
           
            Q1, phi_chain_1, phi_chain_1_numr= single_chain_computation( q_point_1_numr, q_point_5_numr)
            
            ###############################################################################
            ## Chain 2 computation
            
            # generating vtk files to store q_point_2
            vtkfile_phi_chain = File('phi_chain_2.pvd')
            
            Q2, phi_chain_2, phi_chain_2_numr= single_chain_computation( q_point_2_numr, q_point_5_numr)
            
            ###############################################################################
            ## Chain 3 computation
            
            # generating vtk files to store q_point_1
            vtkfile_phi_chain = File('phi_chain_3.pvd')
            
            Q3, phi_chain_3, phi_chain_3_numr= single_chain_computation( q_point_3_numr, q_point_5_numr)        
            
            ###############################################################################
            ## Chain 4 computation
            
            # generating vtk files to store q_point_4
            vtkfile_phi_chain = File('phi_chain_4.pvd')
            
            Q4, phi_chain_4, phi_chain_4_numr= single_chain_computation( q_point_4_numr, q_point_5_numr)   
            
            
            
            
            
            ###############################################################################
            ##getting Q at each chain mid point        
            
            
            
            Q1_mid= Q1[round(n_t/2)] 
            Q2_mid= Q2[round(n_t/2)] 
            Q3_mid= Q3[round(n_t/2)] 
            Q4_mid= Q4[round(n_t/2)] 
            
            
            
            
            #converting phi_MF into fem functions, for w_next
            phi_MF=Function(V)  # phi function
            phi_MF_numr=np.zeros(n_mesh) #phi over mesh nodes(vector of nx.ny.nz size)
            
            phi_MF_numr= phi_chain_1_numr + phi_chain_2_numr + phi_chain_3_numr + phi_chain_4_numr                
            phi_MF.vector()[:]= phi_MF_numr  
            
            #saving phi_MF 
            vtkfile_phi_MF = File('phi_MF.pvd')
            vtkfile_phi_MF << (phi_MF)
            
            
            #total free energy, for w_next
            H_next = k_B*Temp*(1/(2*u0))*assemble(w_next*w_next*dx) - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  )       
            
            H_next_entropic= - k_B*Temp* ( math.log(Q1_mid) + math.log(Q2_mid) + math.log(Q3_mid)+  math.log(Q4_mid)  )
            H_next_interaction=  k_B*Temp*(1/(2*u0))*assemble(w_next*w_next*dx)
            
            # computing relative change in H and H_next
            delta_H_ratio=(H_next-H)/H
            delta_H_ratio_array[count]=delta_H_ratio
            
            print(delta_H_ratio)
            
            if abs((H_next-H)/H) < delta_H_ratio_threshold:
                print('iteration ended after break')
                break
            
            
            ############################################################### 
            ## To get w for next time step
            ###############################################################
            
            
            ############################################################### 
            ## Shifting functions 
            ############################################################### 
            
            phi_MF_B13= interpolate (phi_MF, V_B13)
            File('phi_MF_B13.pvd') << (phi_MF_B13)
        
        
            phi_MF_B1= Function(V_B1)
            phi_MF_B1.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B1.pvd') << (phi_MF_B1)
        
            phi_MF_B2= Function(V_B2)
            phi_MF_B2.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B2.pvd') << (phi_MF_B2)
        
            phi_MF_B3= Function(V_B3)
            phi_MF_B3.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B3.pvd') << (phi_MF_B3)
            
            phi_MF_B4= Function(V_B4)
            phi_MF_B4.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B4.pvd') << (phi_MF_B4)
            
            phi_MF_B5= Function(V_B5)
            phi_MF_B5.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B5.pvd') << (phi_MF_B5)
            
            phi_MF_B6= Function(V_B6)
            phi_MF_B6.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B6.pvd') << (phi_MF_B6)
            
            phi_MF_B7= Function(V_B7)
            phi_MF_B7.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B7.pvd') << (phi_MF_B7)
            
            phi_MF_B8= Function(V_B8)
            phi_MF_B8.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B8.pvd') << (phi_MF_B8)
            
            phi_MF_B9= Function(V_B9)
            phi_MF_B9.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B9.pvd') << (phi_MF_B9)
            
            phi_MF_B10= Function(V_B10)
            phi_MF_B10.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B10.pvd') << (phi_MF_B10)
            
            phi_MF_B11= Function(V_B11)
            phi_MF_B11.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B11.pvd') << (phi_MF_B11)
            
            phi_MF_B12= Function(V_B12)
            phi_MF_B12.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B12.pvd') << (phi_MF_B12)
            
            
            phi_MF_B14= Function(V_B14)
            phi_MF_B14.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B14.pvd') << (phi_MF_B14)
            
            phi_MF_B15= Function(V_B15)
            phi_MF_B15.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B15.pvd') << (phi_MF_B15)
            
            phi_MF_B16= Function(V_B16)
            phi_MF_B16.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B16.pvd') << (phi_MF_B16)
            
            phi_MF_B17= Function(V_B17)
            phi_MF_B17.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B17.pvd') << (phi_MF_B17)
            
            phi_MF_B18= Function(V_B18)
            phi_MF_B18.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B18.pvd') << (phi_MF_B18)
            
            phi_MF_B19= Function(V_B19)
            phi_MF_B19.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B19.pvd') << (phi_MF_B19)
            
            phi_MF_B20= Function(V_B20)
            phi_MF_B20.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B20.pvd') << (phi_MF_B20)
            
            phi_MF_B21= Function(V_B21)
            phi_MF_B21.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B21.pvd') << (phi_MF_B21)
            
            phi_MF_B22= Function(V_B22)
            phi_MF_B22.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B22.pvd') << (phi_MF_B22)
            
            phi_MF_B23= Function(V_B23)
            phi_MF_B23.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B23.pvd') << (phi_MF_B23)
            
            phi_MF_B24= Function(V_B24)
            phi_MF_B24.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B24.pvd') << (phi_MF_B24)
            
            phi_MF_B25= Function(V_B25)
            phi_MF_B25.vector()[:]= phi_MF_B13.vector()[:]
            File('phi_MF_B25.pvd') << (phi_MF_B25)
        
        

        
            print('stage 3 done')
        
   
        
        
            ############################################################### 
            ## Getting solutions on submesh 
            
            phi_MF_B13_submesh= interpolate (phi_MF_B13, V_B13_submesh)
            File('phi_MF_B13_submesh.pvd') << (phi_MF_B13_submesh)
        
        
            phi_MF_B1_submesh= interpolate (phi_MF_B1, V_B1_submesh)
            File('phi_MF_B1_submesh.pvd') << (phi_MF_B1_submesh)
        
        
            phi_MF_B2_submesh= interpolate (phi_MF_B2, V_B2_submesh)
            File('phi_MF_B2_submesh.pvd') << (phi_MF_B2_submesh)
            
            
            phi_MF_B3_submesh= interpolate (phi_MF_B3, V_B3_submesh)
            File('phi_MF_B3_submesh.pvd') << (phi_MF_B3_submesh)
            
            
            phi_MF_B4_submesh= interpolate (phi_MF_B4, V_B4_submesh)
            File('phi_MF_B4_submesh.pvd') << (phi_MF_B4_submesh)
            
            
            phi_MF_B5_submesh= interpolate (phi_MF_B5, V_B5_submesh)
            File('phi_MF_B5_submesh.pvd') << (phi_MF_B5_submesh)
            
            
            phi_MF_B6_submesh= interpolate (phi_MF_B6, V_B6_submesh)
            File('phi_MF_B6_submesh.pvd') << (phi_MF_B6_submesh)
            
            
            phi_MF_B7_submesh= interpolate (phi_MF_B7, V_B7_submesh)
            File('phi_MF_B7_submesh.pvd') << (phi_MF_B7_submesh)
            
            
            phi_MF_B8_submesh= interpolate (phi_MF_B8, V_B8_submesh)
            File('phi_MF_B8_submesh.pvd') << (phi_MF_B8_submesh)
            
            
            phi_MF_B9_submesh= interpolate (phi_MF_B9, V_B9_submesh)
            File('phi_MF_B9_submesh.pvd') << (phi_MF_B9_submesh)
            
            
            phi_MF_B10_submesh= interpolate (phi_MF_B10, V_B10_submesh)
            File('phi_MF_B10_submesh.pvd') << (phi_MF_B10_submesh)
            
            
            phi_MF_B11_submesh= interpolate (phi_MF_B11, V_B11_submesh)
            File('phi_MF_B11_submesh.pvd') << (phi_MF_B11_submesh)
            
            
            phi_MF_B12_submesh= interpolate (phi_MF_B12, V_B12_submesh)
            File('phi_MF_B12_submesh.pvd') << (phi_MF_B12_submesh)
            
            
            
            phi_MF_B14_submesh= interpolate (phi_MF_B14, V_B14_submesh)
            File('phi_MF_B14_submesh.pvd') << (phi_MF_B14_submesh)
            
            
            phi_MF_B15_submesh= interpolate (phi_MF_B15, V_B15_submesh)
            File('phi_MF_B15_submesh.pvd') << (phi_MF_B15_submesh)
            
            
            phi_MF_B16_submesh= interpolate (phi_MF_B16, V_B16_submesh)
            File('phi_MF_B16_submesh.pvd') << (phi_MF_B16_submesh)
            
            
            phi_MF_B17_submesh= interpolate (phi_MF_B17, V_B17_submesh)
            File('phi_MF_B17_submesh.pvd') << (phi_MF_B17_submesh)
            
            
            phi_MF_B18_submesh= interpolate (phi_MF_B18, V_B18_submesh)
            File('phi_MF_B18_submesh.pvd') << (phi_MF_B18_submesh)
            
            
            phi_MF_B19_submesh= interpolate (phi_MF_B19, V_B19_submesh)
            File('phi_MF_B19_submesh.pvd') << (phi_MF_B19_submesh)
            
            
            phi_MF_B20_submesh= interpolate (phi_MF_B20, V_B20_submesh)
            File('phi_MF_B20_submesh.pvd') << (phi_MF_B20_submesh)
            
            
            phi_MF_B21_submesh= interpolate (phi_MF_B21, V_B21_submesh)
            File('phi_MF_B21_submesh.pvd') << (phi_MF_B21_submesh)
            
            
            phi_MF_B22_submesh= interpolate (phi_MF_B22, V_B22_submesh)
            File('phi_MF_B22_submesh.pvd') << (phi_MF_B22_submesh)
            
            
            phi_MF_B23_submesh= interpolate (phi_MF_B23, V_B23_submesh)
            File('phi_MF_B23_submesh.pvd') << (phi_MF_B23_submesh)
            
            
            phi_MF_B24_submesh= interpolate (phi_MF_B24, V_B24_submesh)
            File('phi_MF_B24_submesh.pvd') << (phi_MF_B24_submesh)
            
            
            phi_MF_B25_submesh= interpolate (phi_MF_B25, V_B25_submesh)
            File('phi_MF_B25_submesh.pvd') << (phi_MF_B25_submesh)
        
        
        
        
        
            print('stage 4 done')
        
        
        
        
        
            ####################################### converting submesh functions to big-box-mesh function    
            
            
            ##################### getting phi_MF_B13_submesh on V_big
        
            dof_v_V_B13_submesh= np.array(dof_to_vertex_map(V_B13_submesh), dtype=int) 
            v_V_B13_submesh = mesh_B13_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B13_submesh_on_V_big = Function(V_big)
            array = phi_MF_B13_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B13_submesh[dof_v_V_B13_submesh]]] = phi_MF_B13_submesh.vector()[:]
            phi_MF_B13_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B13_submesh_on_V_big.pvd') << (phi_MF_B13_submesh_on_V_big)
        
            
            ##################### getting phi_MF_B1_submesh on V_big
        
            dof_v_V_B1_submesh= np.array(dof_to_vertex_map(V_B1_submesh), dtype=int) 
            v_V_B1_submesh = mesh_B1_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B1_submesh_on_V_big = Function(V_big)
            array = phi_MF_B1_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B1_submesh[dof_v_V_B1_submesh]]] = phi_MF_B1_submesh.vector()[:]
            phi_MF_B1_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B1_submesh_on_V_big.pvd') << (phi_MF_B1_submesh_on_V_big)
            
            ##################### getting phi_MF_B2_submesh on V_big
        
            dof_v_V_B2_submesh = np.array(dof_to_vertex_map(V_B2_submesh), dtype=int) 
            v_V_B2_submesh = mesh_B2_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B2_submesh_on_V_big = Function(V_big)
            array = phi_MF_B2_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B2_submesh[dof_v_V_B2_submesh]]] = phi_MF_B2_submesh.vector()[:]
            phi_MF_B2_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B2_submesh_on_V_big.pvd') << (phi_MF_B2_submesh_on_V_big)
            
            ##################### getting phi_MF_B3_submesh on V_big
        
            dof_v_V_B3_submesh= np.array(dof_to_vertex_map(V_B3_submesh), dtype=int) 
            v_V_B3_submesh = mesh_B3_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B3_submesh_on_V_big = Function(V_big)
            array = phi_MF_B3_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B3_submesh[dof_v_V_B3_submesh]]] = phi_MF_B3_submesh.vector()[:]
            phi_MF_B3_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B3_submesh_on_V_big.pvd') << (phi_MF_B3_submesh_on_V_big)
            
            ##################### getting phi_MF_B4_submesh on V_big
        
            dof_v_V_B4_submesh= np.array(dof_to_vertex_map(V_B4_submesh), dtype=int) 
            v_V_B4_submesh = mesh_B4_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B4_submesh_on_V_big = Function(V_big)
            array = phi_MF_B4_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B4_submesh[dof_v_V_B4_submesh]]] = phi_MF_B4_submesh.vector()[:]
            phi_MF_B4_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B4_submesh_on_V_big.pvd') << (phi_MF_B4_submesh_on_V_big)
            
            ##################### getting phi_MF_B5_submesh on V_big
        
            dof_v_V_B5_submesh= np.array(dof_to_vertex_map(V_B5_submesh), dtype=int) 
            v_V_B5_submesh = mesh_B5_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B5_submesh_on_V_big = Function(V_big)
            array = phi_MF_B5_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B5_submesh[dof_v_V_B5_submesh]]] = phi_MF_B5_submesh.vector()[:]
            phi_MF_B5_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B5_submesh_on_V_big.pvd') << (phi_MF_B5_submesh_on_V_big)
            
            ##################### getting phi_MF_B6_submesh on V_big
        
            dof_v_V_B6_submesh= np.array(dof_to_vertex_map(V_B6_submesh), dtype=int) 
            v_V_B6_submesh = mesh_B6_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B6_submesh_on_V_big = Function(V_big)
            array = phi_MF_B6_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B6_submesh[dof_v_V_B6_submesh]]] = phi_MF_B6_submesh.vector()[:]
            phi_MF_B6_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B6_submesh_on_V_big.pvd') << (phi_MF_B6_submesh_on_V_big)
            
            ##################### getting phi_MF_B7_submesh on V_big
        
            dof_v_V_B7_submesh= np.array(dof_to_vertex_map(V_B7_submesh), dtype=int) 
            v_V_B7_submesh = mesh_B7_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B7_submesh_on_V_big = Function(V_big)
            array = phi_MF_B7_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B7_submesh[dof_v_V_B7_submesh]]] = phi_MF_B7_submesh.vector()[:]
            phi_MF_B7_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B7_submesh_on_V_big.pvd') << (phi_MF_B7_submesh_on_V_big)
            
            ##################### getting phi_MF_B8_submesh on V_big
        
            dof_v_V_B8_submesh= np.array(dof_to_vertex_map(V_B8_submesh), dtype=int) 
            v_V_B8_submesh = mesh_B8_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B8_submesh_on_V_big = Function(V_big)
            array = phi_MF_B8_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B8_submesh[dof_v_V_B8_submesh]]] = phi_MF_B8_submesh.vector()[:]
            phi_MF_B8_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B8_submesh_on_V_big.pvd') << (phi_MF_B8_submesh_on_V_big)
            
            ##################### getting phi_MF_B9_submesh on V_big
        
            dof_v_V_B9_submesh= np.array(dof_to_vertex_map(V_B9_submesh), dtype=int) 
            v_V_B9_submesh = mesh_B9_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B9_submesh_on_V_big = Function(V_big)
            array = phi_MF_B9_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B9_submesh[dof_v_V_B9_submesh]]] = phi_MF_B9_submesh.vector()[:]
            phi_MF_B9_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B9_submesh_on_V_big.pvd') << (phi_MF_B9_submesh_on_V_big)
            
            ##################### getting phi_MF_B10_submesh on V_big
        
            dof_v_V_B10_submesh= np.array(dof_to_vertex_map(V_B10_submesh), dtype=int) 
            v_V_B10_submesh = mesh_B10_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B10_submesh_on_V_big = Function(V_big)
            array = phi_MF_B10_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B10_submesh[dof_v_V_B10_submesh]]] = phi_MF_B10_submesh.vector()[:]
            phi_MF_B10_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B10_submesh_on_V_big.pvd') << (phi_MF_B10_submesh_on_V_big)
            
            ##################### getting phi_MF_B11_submesh on V_big
        
            dof_v_V_B11_submesh= np.array(dof_to_vertex_map(V_B11_submesh), dtype=int) 
            v_V_B11_submesh = mesh_B11_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B11_submesh_on_V_big = Function(V_big)
            array = phi_MF_B11_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B11_submesh[dof_v_V_B11_submesh]]] = phi_MF_B11_submesh.vector()[:]
            phi_MF_B11_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B11_submesh_on_V_big.pvd') << (phi_MF_B11_submesh_on_V_big)
            
            ##################### getting phi_MF_B12_submesh on V_big
        
            dof_v_V_B12_submesh= np.array(dof_to_vertex_map(V_B12_submesh), dtype=int) 
            v_V_B12_submesh = mesh_B12_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B12_submesh_on_V_big = Function(V_big)
            array = phi_MF_B12_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B12_submesh[dof_v_V_B12_submesh]]] = phi_MF_B12_submesh.vector()[:]
            phi_MF_B12_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B12_submesh_on_V_big.pvd') << (phi_MF_B12_submesh_on_V_big)
            
            ##################### getting phi_MF_B13_submesh on V_big
        
            dof_v_V_B13_submesh= np.array(dof_to_vertex_map(V_B13_submesh), dtype=int) 
            v_V_B13_submesh = mesh_B13_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B13_submesh_on_V_big = Function(V_big)
            array = phi_MF_B13_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B13_submesh[dof_v_V_B13_submesh]]] = phi_MF_B13_submesh.vector()[:]
            phi_MF_B13_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B13_submesh_on_V_big.pvd') << (phi_MF_B13_submesh_on_V_big)
            
            ##################### getting phi_MF_B14_submesh on V_big
        
            dof_v_V_B14_submesh= np.array(dof_to_vertex_map(V_B14_submesh), dtype=int) 
            v_V_B14_submesh = mesh_B14_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B14_submesh_on_V_big = Function(V_big)
            array = phi_MF_B14_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B14_submesh[dof_v_V_B14_submesh]]] = phi_MF_B14_submesh.vector()[:]
            phi_MF_B14_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B14_submesh_on_V_big.pvd') << (phi_MF_B14_submesh_on_V_big)
            
            ##################### getting phi_MF_B15_submesh on V_big
        
            dof_v_V_B15_submesh= np.array(dof_to_vertex_map(V_B15_submesh), dtype=int) 
            v_V_B15_submesh = mesh_B15_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B15_submesh_on_V_big = Function(V_big)
            array = phi_MF_B15_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B15_submesh[dof_v_V_B15_submesh]]] = phi_MF_B15_submesh.vector()[:]
            phi_MF_B15_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B15_submesh_on_V_big.pvd') << (phi_MF_B15_submesh_on_V_big)
            
            ##################### getting phi_MF_B16_submesh on V_big
        
            dof_v_V_B16_submesh= np.array(dof_to_vertex_map(V_B16_submesh), dtype=int) 
            v_V_B16_submesh = mesh_B16_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B16_submesh_on_V_big = Function(V_big)
            array = phi_MF_B16_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B16_submesh[dof_v_V_B16_submesh]]] = phi_MF_B16_submesh.vector()[:]
            phi_MF_B16_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B16_submesh_on_V_big.pvd') << (phi_MF_B16_submesh_on_V_big)
            
            ##################### getting phi_MF_B17_submesh on V_big
        
            dof_v_V_B17_submesh= np.array(dof_to_vertex_map(V_B17_submesh), dtype=int) 
            v_V_B17_submesh = mesh_B17_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B17_submesh_on_V_big = Function(V_big)
            array = phi_MF_B17_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B17_submesh[dof_v_V_B17_submesh]]] = phi_MF_B17_submesh.vector()[:]
            phi_MF_B17_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B17_submesh_on_V_big.pvd') << (phi_MF_B17_submesh_on_V_big)
            
            ##################### getting phi_MF_B18_submesh on V_big
        
            dof_v_V_B18_submesh= np.array(dof_to_vertex_map(V_B18_submesh), dtype=int) 
            v_V_B18_submesh = mesh_B18_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B18_submesh_on_V_big = Function(V_big)
            array = phi_MF_B18_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B18_submesh[dof_v_V_B18_submesh]]] = phi_MF_B18_submesh.vector()[:]
            phi_MF_B18_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B18_submesh_on_V_big.pvd') << (phi_MF_B18_submesh_on_V_big)
            
            ##################### getting phi_MF_B19_submesh on V_big
        
            dof_v_V_B19_submesh= np.array(dof_to_vertex_map(V_B19_submesh), dtype=int) 
            v_V_B19_submesh = mesh_B19_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B19_submesh_on_V_big = Function(V_big)
            array = phi_MF_B19_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B19_submesh[dof_v_V_B19_submesh]]] = phi_MF_B19_submesh.vector()[:]
            phi_MF_B19_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B19_submesh_on_V_big.pvd') << (phi_MF_B19_submesh_on_V_big)
            
            ##################### getting phi_MF_B20_submesh on V_big
        
            dof_v_V_B20_submesh= np.array(dof_to_vertex_map(V_B20_submesh), dtype=int) 
            v_V_B20_submesh = mesh_B20_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B20_submesh_on_V_big = Function(V_big)
            array = phi_MF_B20_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B20_submesh[dof_v_V_B20_submesh]]] = phi_MF_B20_submesh.vector()[:]
            phi_MF_B20_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B20_submesh_on_V_big.pvd') << (phi_MF_B20_submesh_on_V_big)
            
            ##################### getting phi_MF_B21_submesh on V_big
        
            dof_v_V_B21_submesh= np.array(dof_to_vertex_map(V_B21_submesh), dtype=int) 
            v_V_B21_submesh = mesh_B21_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B21_submesh_on_V_big = Function(V_big)
            array = phi_MF_B21_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B21_submesh[dof_v_V_B21_submesh]]] = phi_MF_B21_submesh.vector()[:]
            phi_MF_B21_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B21_submesh_on_V_big.pvd') << (phi_MF_B21_submesh_on_V_big)
            
            ##################### getting phi_MF_B22_submesh on V_big
        
            dof_v_V_B22_submesh= np.array(dof_to_vertex_map(V_B22_submesh), dtype=int) 
            v_V_B22_submesh = mesh_B22_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B22_submesh_on_V_big = Function(V_big)
            array = phi_MF_B22_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B22_submesh[dof_v_V_B22_submesh]]] = phi_MF_B22_submesh.vector()[:]
            phi_MF_B22_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B22_submesh_on_V_big.pvd') << (phi_MF_B22_submesh_on_V_big)
            
            ##################### getting phi_MF_B23_submesh on V_big
        
            dof_v_V_B23_submesh= np.array(dof_to_vertex_map(V_B23_submesh), dtype=int) 
            v_V_B23_submesh = mesh_B23_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B23_submesh_on_V_big = Function(V_big)
            array = phi_MF_B23_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B23_submesh[dof_v_V_B23_submesh]]] = phi_MF_B23_submesh.vector()[:]
            phi_MF_B23_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B23_submesh_on_V_big.pvd') << (phi_MF_B23_submesh_on_V_big)
            
            ##################### getting phi_MF_B24_submesh on V_big
        
            dof_v_V_B24_submesh= np.array(dof_to_vertex_map(V_B24_submesh), dtype=int) 
            v_V_B24_submesh = mesh_B24_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B24_submesh_on_V_big = Function(V_big)
            array = phi_MF_B24_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B24_submesh[dof_v_V_B24_submesh]]] = phi_MF_B24_submesh.vector()[:]
            phi_MF_B24_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B24_submesh_on_V_big.pvd') << (phi_MF_B24_submesh_on_V_big)
            
            ##################### getting phi_MF_B25_submesh on V_big
        
            dof_v_V_B25_submesh= np.array(dof_to_vertex_map(V_B25_submesh), dtype=int) 
            v_V_B25_submesh = mesh_B25_submesh.data().array('parent_vertex_indices', 0)
            
            v_dof_V_big = np.array(vertex_to_dof_map(V_big), dtype=int)
            
            phi_MF_B25_submesh_on_V_big = Function(V_big)
            array = phi_MF_B25_submesh_on_V_big.vector()[:]
            # Compose the maps dof of Vr --> vertices of Vr --> vertices of V --> dof of V
            array[v_dof_V_big[v_V_B25_submesh[dof_v_V_B25_submesh]]] = phi_MF_B25_submesh.vector()[:]
            phi_MF_B25_submesh_on_V_big.vector()[:] = array
            
            File('phi_MF_B25_submesh_on_V_big.pvd') << (phi_MF_B25_submesh_on_V_big)
        
        
        
        
        
            print('stage 5 done')
        
        
        
            
            
            ############################################################### 
            ## getting w for next time step
            ############################################################### 
            
            phi_MF_on_V_big= Function(V_big)
            phi_MF_on_V_big.vector()[:] = phi_MF_B1_submesh_on_V_big.vector()[:] +  phi_MF_B2_submesh_on_V_big.vector()[:] +  phi_MF_B3_submesh_on_V_big.vector()[:] +  phi_MF_B4_submesh_on_V_big.vector()[:] +  phi_MF_B5_submesh_on_V_big.vector()[:] + phi_MF_B6_submesh_on_V_big.vector()[:] +  phi_MF_B7_submesh_on_V_big.vector()[:] + phi_MF_B8_submesh_on_V_big.vector()[:] +  phi_MF_B9_submesh_on_V_big.vector()[:]  + phi_MF_B10_submesh_on_V_big.vector()[:] + phi_MF_B11_submesh_on_V_big.vector()[:] + phi_MF_B12_submesh_on_V_big.vector()[:] + phi_MF_B13_submesh_on_V_big.vector()[:] + phi_MF_B14_submesh_on_V_big.vector()[:] + phi_MF_B15_submesh_on_V_big.vector()[:] + phi_MF_B16_submesh_on_V_big.vector()[:] + phi_MF_B17_submesh_on_V_big.vector()[:] + phi_MF_B18_submesh_on_V_big.vector()[:] + phi_MF_B19_submesh_on_V_big.vector()[:] + phi_MF_B20_submesh_on_V_big.vector()[:] + phi_MF_B21_submesh_on_V_big.vector()[:] + phi_MF_B22_submesh_on_V_big.vector()[:] + phi_MF_B23_submesh_on_V_big.vector()[:] + phi_MF_B24_submesh_on_V_big.vector()[:] + phi_MF_B25_submesh_on_V_big.vector()[:]  
        
            
            #saving phi_MF on V_big
            vtkfile_phi_MF = File('phi_MF_on_V_big.pvd')
            vtkfile_phi_MF << (phi_MF_on_V_big)
            
            ########## getting phi_MF on V
            phi_MF_on_V_big_to_V= interpolate (phi_MF_on_V_big, V)
            File('phi_MF_on_V_big_to_V.pvd') << (phi_MF_on_V_big_to_V)
            
            
            #defining w field for next step
            w=Function(V)  # phi function
            
            w.vector()[:]= u0* phi_MF_on_V_big_to_V.vector()[:] 
            # w=w_next
    
        ############################ End of iteration for finding mean field w
        
        # saving free energy values 
        W[i1, i2]=H_next #total free energy
        W_entropic[i1, i2]=H_next_entropic #entropic free energy
        W_interaction[i1, i2]=H_next_interaction #interaction free energy


############################ End of computation





############################ 
## Plotting
############################ 


# ## checking plot for q initial condition
# fig = plt.figure()                                                               
# ax = fig.add_subplot(111, projection='3d')                                       
# ax.scatter(dof_x, dof_y, q_point_7_numr[0][:], c='b', marker='.')                  
# plt.show() 

# ## checking plot for q_star initial condition
# fig = plt.figure()                                                               
# ax = fig.add_subplot(111, projection='3d')                                       
# ax.scatter(dof_x, dof_y, q_point_4_numr[0][:], c='b', marker='.')                  
# plt.show() 

# ## checking plot for w 
# fig = plt.figure()                                                               
# ax = fig.add_subplot(111, projection='3d')                                       
# ax.scatter(dof_x, dof_y, w_next.vector().get_local(), c='b', marker='.')                  
# plt.show()

# # checking plot for phi_MF
# fig = plt.figure()                                                               
# ax = fig.add_subplot(111, projection='3d')                                       
# ax.scatter(dof_x, dof_y, phi_MF.vector().get_local(), c='b', marker='.')                  
# plt.show() 










