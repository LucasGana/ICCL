# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:56:15 2023

@author: lucas
"""
#%%
import numpy as np
import gurobipy as gp
from gurobipy import GRB
#%%
#PARAMETERS
# total number of periods in the horizon
period = [0,1,2,3,4]

#Number of item
item = [0,1,2]
    
#Flow returns in period t for item i
r  = [[0,180,180,180,0],
      [0,180,180,180,0],
      [0,180,180,180,0]]

#Demand in period t for item t
d = [[0,20,20,20,20],
     [0,20,20,20,20],
     [0,20,20,20,20]]


#Lost sales cost
p = 200

#Production Ratio
prodRatio = 2

#Setup Ratio
setupRatio = 1000

#Seed
np.random.seed(345)

# unit inventory serviceable holding cost in period t.This is charged for inventory at the end of period (h_t)
hs_t = np.random.uniform(7, 12, (len(item), len(period)))
# unit production (procurement) cost (variable cost) in period t (c_t)
########PARA LA VERSIÓN JOURNAL AGREGAR H BARRITA######
c_t = np.random.uniform(0.8*prodRatio, 1.2*prodRatio, (len(item), len(period)))
# setup cost (fixed cost) of production(procurement) in period t (s_t)
s_t = np.random.uniform(0.8*setupRatio, 1.2*setupRatio, (len(item), len(period)))
# Unit inventory remanufacturing holding cost in period t.
hr_t = np.random.uniform(2, 7, (len(item), len(period)))
#Changeovers cost
change = np.random.uniform(0.8*setupRatio, 1.2*setupRatio, (len(item), len(item), len(period)))
#Alpha
alpha = np.random.uniform(1, 6, (len(item), len(item)))

#Phi
phi = {}
for i in item:
    for k in range(len(period)):
        phi[k, i] = sum(d[i][v] for v in range(k, len(period)))
        for j in item:
            if phi[k, i] > sum(r[j][l]/alpha[i][j] for l in range(k+1)):
                phi[k, i] = sum(r[j][l]/alpha[i][j] for l in range(k+1))

#%%
# Model
m = gp.Model("LSP")
m.setParam(GRB.Param.TimeLimit, 300)
# m.setParam(GRB.Param.Cuts,-1)
# VARS
# Quantity production variable
X = m.addVars(
    period, item, lb=0, vtype=GRB.CONTINUOUS, name="Production Quantity")

# Binary Variable in period t
Y = m.addVars(period, item, lb=0,
              vtype=GRB.BINARY, name="Decision")

# Lost Sales Variables in period t
L = m.addVars(period, item, lb=0,
              vtype=GRB.CONTINUOUS, name="Lost Sales")

# Remanufacturing Inventory level variable in period t
SR = m.addVars(
    period, item, lb=0, vtype=GRB.CONTINUOUS, name="Returns Inventory")

# Serviceable Inventory level Variable in period t
SS = m.addVars(
    period, item, lb=0, vtype=GRB.CONTINUOUS, name="Serviceables Inventory")

# Changeovers variables from item i to j in period t
W = m.addVars(period, item, item,
              vtype=GRB.BINARY, name="Changeover")

# %%CONSTRAINTS
# Returns Inventory balance flow constraint to first period
m.addConstrs((SR[t, j] == SR[t-1, j] + r[j][t] - sum(alpha[i][j]*X[t, i] for i in item)
             for t in range(1, len(period)) for j in item), name="Returns Inventory Balance Flow")
m.addConstrs((SR[t, j] ==  r[j][t] - sum(alpha[i][j]*X[t, i] for i in item)
             for t in range(0, 1) for j in item), name="Returns Inventory Balance Flow")

# Inventory balance flow constraint to first period
m.addConstrs((SS[t, i] == SS[t-1, i] + X[t, i] - d[i][t] + L[t, i] for t in range(
    1, len(period)) for i in item), name="Inventory Balance Flow")
m.addConstrs((SS[t, i] == X[t, i] - d[i][t] + L[t, i]
             for t in range(0, 1) for i in item), name="Inventory Balance Flow")

# If decided producted, Y = 1.
m.addConstrs((X[t, i] <= phi[t,i]*Y[t, i]
                  for t in period for i in item), name="Logical Constraints")

# The Lost sales must be under than demand in period t
m.addConstrs(
    (L[t, i] <= d[i][t] for t in period for i in item), name="Lost Sales Constraints")

# the changeovers constraints
m.addConstrs((W[t, i, j] >= Y[t-1, i] + Y[t, j] - 1 for i in item for j in item for t in range(
    1, len(period)) if i != j), name="changeovers constraints")

# just one of the items must be produced
m.addConstrs((sum(Y[t, i] for i in item) <= 1 for t in range(
    len(period))), name="Just one")

# %%OBJECTIVE FUNCTION

m.setObjective((sum(X[t, i]*c_t[i][t] + Y[t, i]*s_t[i][t] + SR[t, i]*hr_t[i][t] + SS[t, i]*hs_t[i][t] + L[t, i]*p for t in range(len(period)) for i in item) +
                sum(W[t, i, j] * change[i, j, t] for (t, i, j) in W.keys() if i != j)), GRB.MINIMIZE)

#Resolve the model
m.optimize()
