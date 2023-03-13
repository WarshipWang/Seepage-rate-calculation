# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 09:16:15 2023
This code is to solve the speeage rate of Hugo A. Loaiciga's model'
Ref: Seepage Face in Steady-State Groundwater Flow Between Two Water Bodies
@author: JuWang
"""
from scipy.optimize import fsolve
import math

#### -----Model Inputs ----


h2            =  29                              # Downstream water head in height, m
h1            =  30                              # Upstream water head in height, m                        
H             =  32                              # Crest height, m
alpha_degree  =  26.5                            # Angle of seepage phase, degree
w             =  500                             # Crest Top Width, m
K             =  1                               # Hydraulic conductivity of Crest, m/d
alpha_radians =  math.radians(alpha_degree)      
tan_theta     =  math.tan(alpha_radians)                             
Cot           =  1/tan_theta
   

###  initial_guess ------
hD_i   = h2  +0.95*(h1-h2)
a_i    = 0.003


###  model solve   -------
def equations(vars):
    hD, a = vars
    L = w + Cot *(H - (h2 + a))
    eq1 = hD - (((h2+a)**2 + (2*a*L)/Cot * (1 + math.log((h2+a)/a)))**0.5)
    eq2 = ((h1 - hD)/Cot * math.log(H/(H-hD)) - a/Cot * (1 + math.log((h2+a)/a)))
    return [eq1, eq2]

hD, a = fsolve(equations, (hD_i, a_i))
L = w + Cot *(H - (h2 + a))


q = K * (hD**2 - (h2 + a)**2) / (2 * L)        # outflow rate from seepage phase, m3/m/d    
 
print(f"hD = {hD:.2f}, a = {a:.3f}, q = {q:.2f}")
