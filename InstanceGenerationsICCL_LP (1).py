# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 20:14:29 2023

@author: lucas
"""

import numpy as np
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import time

########FUNCTIONS#######
# %%WagnerWithinInequalities
def WagnerWithin(model, where):
    if where == GRB.Callback.MIPNODE:
        t1 = time.time()
        obj_bound = model.cbGet(gp.GRB.Callback.MIPNODE_NODCNT)
        if obj_bound == 0:
            m._ContCallback += 1
            m._BestBound = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
            # vals_y = model.cbGetNodeRel(Y)
            # vals_ss = model.cbGetNodeRel(SS)
            # vals_l = model.cbGetNodeRel(L)

            # for i in item:
            #     for k in period:
            #         for l in range(k, len(period)):
            #             if k == 0:
            #                 left = sum(
            #                     vals_y[t, i] * sum(d[i][v] for v in range(t, l+1)) for t in range(k, l+1))
            #                 right = sum(d[i][t] - vals_l[t, i]
            #                             for t in range(k, l+1))
            #                 if left < right:
            #                     m._ContCut += 1
            #                     model.cbCut(sum(Y[t, i] * sum(d[i][v] for v in range(t, l+1)) for t in range(k, l+1))
            #                                 >= sum(d[i][t] - L[t, i] for t in range(k, l+1)))
            #             else:
            #                 left = vals_ss[k-1, i] + sum(vals_y[t, i] * sum(
            #                     d[i][v] for v in range(t, l+1)) for t in range(k, l+1))
            #                 right = sum(d[i][t] - vals_l[t, i]
            #                             for t in range(k, l+1))
            #                 if left < right:
            #                     m._ContCut += 1
            #                     model.cbCut(SS[k-1, i] + sum(Y[t, i] * sum(d[i][v] for v in range(t, l+1)) for t in range(k, l+1))
            #                                 >= sum(d[i][t] - L[t, i] for t in range(k, l+1)))
        else:
            
            m._time += time.time() - t1
            m.terminate()
# %% SINGLE ITEM
def SingleItem(m):
    U = np.zeros((len(period), len(period), len(item)))
    for i in item:
        for k in period:
            for t in range(k, len(period)):
                left = d[i][t]*(1-sum(Y[v, i].X
                                for v in range(k, t+1))) - L[t, i].X
                if left > 0:
                    U[k, t, i] = 1
    t_max = {}                    
    for i in item:
        for k in period[0:len(period)]:
            t_max[k, i] = k-1
            for t in range(k, len(period)):
                if U[k, t, i] == 1:
                    t_max[k, i] = t

    for i in item:
        for k in period[0:len(period)]:
            if k == 0:

                left = sum(Y[t, i].X * sum(d[i][v]*U[k][v][i]
                           for v in range(t, len(period))) for t in range(k, len(period)))
                right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                            for t in range(k, len(period)))
                if left + 0.01 < right:
                    m._ContCut += 1

                    m.addConstr(sum(Y[t, i] * sum(d[i][v]*U[k][v][i] for v in range(t, len(period))) for t in range(k, t_max[k, i]+1))
                                >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, len(period))), name="Cuts")
                   
                    
            else:

                left = SS[k-1, i].X + sum(Y[t, i].X * sum(
                    d[i][v]*U[k][v][i] for v in range(t, len(period))) for t in range(k, t_max[k, i]+1))
                right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                            for t in range(k, len(period)))
                if left + 0.01 < right:
                    m._ContCut += 1
                    m.addConstr(SS[k-1, i] + sum(Y[t, i] * sum(d[i][v]*U[k][v][i] for v in range(t, len(period))) for t in range(k, len(period)))
                                >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, len(period))), name="Cuts")
                    
                    
# %% (k,U) valid inequalities
def SingleItemAlpha(model):
    U = np.zeros((len(period), len(period), len(item)))
    for i in item:
        for k in period:
            for t in range(k, len(period)):
                left = d[i][t]*(1-sum(Y[v, i].X
                                for v in range(k, t+1))) - L[t, i].X
                if left > 0:
                    U[k, t, i] = 1
            
    phi = np.zeros((len(period), len(item)))
    t_max = {}
    for i in item:
        for k in period[0:len(period)]:
            t_max[k, i] = k-1
            for t in range(k, len(period)):
                if U[k, t, i] == 1:
                    t_max[k, i] = t

    for i in item:
        for k in period[0:len(period)]:
            min_alpha = np.zeros(len(item))
            phi[k, i] = 100000000
            for j in item:
                min_alpha[j] = sum(r[j][l]/alpha[i][j]
                                   for l in range(k+1))
                if phi[k, i] > min_alpha[j]:
                    phi[k, i] = min_alpha[j]
            
    for i in item:
        for k in period[0:len(period)]:
            if k == 0:

                left = sum(Y[t, i].X * min(phi[t, i], sum(d[i][v]*U[k][v][i]
                           for v in range(t, len(period)))) for t in range(k, t_max[k, i]+1))
                right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                            for t in range(k, len(period)))
                if left + 0.01 < right:
                    m._ContCut += 1
                    model.addConstr(sum(Y[t, i] * min(phi[t, i], sum(d[i][v]*U[k][v][i] for v in range(t, len(period)))) for t in range(k, t_max[k, i]+1))
                                >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, len(period))))
                    
            else:
                left = SS[k-1, i].X + sum(Y[t, i].X * min(phi[t, i], sum(d[i][v]*U[k][v][i] for v in range(t, len(period)))) for t in range(k, t_max[k, i]+1))
                
                right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                            for t in range(k, len(period)))
                if left + 0.01 < right:
                    m._ContCut += 1
                    model.addConstr(SS[k-1, i] + sum(Y[t, i] * min(phi[t, i], sum(d[i][v]*U[k][v][i] for v in range(t, len(period)))) for t in range(k, t_max[k, i]+1))
                                >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, len(period))))
                    
                        
        
# %% (l,k,U) valid inequalities
def SingleItemReturnsa(model):
    m._ContCallback += 1
     
    U = np.zeros((len(period), len(period), len(item)))
    for i in item:
        for k in period:
            for t in range(k, len(period)):
                left = d[i][t]*(1-sum(Y[v, i].X
                                for v in range(k, t+1))) - L[t, i].X
                if left > 0:
                    U[k, t, i] = 1
   
    t_max = {}
    for i in item:
        for k in period[0:len(period)]:
            t_max[k, i] = k-1
            for t in range(k, len(period)):
                if U[k, t, i] == 1:
                    t_max[k, i] = t

    for i in item:
        for j in item:
            for k in period[0:len(period)]:
                for l in range(0, k+1):
                    if l == 0 and k == 0:
                        left = sum(Y[t, i].X * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                   for t in range(k, t_max[k, i]+1))
                        right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                                    for t in range(k, t_max[k, i]+1))

                        if left + 0.01 < right:
                            m._ContCut += 1
                            model.addConstr(sum(Y[t, i] * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                            for t in range(k, t_max[k, i]+1))
                                        >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, t_max[k, i]+1)))
                    if l > 0 and k == 0:
                        left = SR[l-1, j].X/alpha[i][j] + sum(Y[t, i].X * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                                     for t in range(k, t_max[k, i]+1))
                        right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                                    for t in range(k, t_max[k, i]+1))

                        if left + 0.01 < right:
                            m._ContCut += 1
                            model.addConstr(SR[l-1, j]/alpha[i][j] + sum(Y[t, i] * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                                         for t in range(k, t_max[k, i]+1))
                                        >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, t_max[k, i]+1)))
                    if l == 0 and k > 0:
                        left = SS[k-1, i].X + sum(Y[t, i].X * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                                     for t in range(k, t_max[k, i]+1))
                        right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                                    for t in range(k, t_max[k, i]+1))

                        if left + 0.01 < right:
                            m._ContCut += 1
                            model.addConstr(SS[k-1, i] + sum(Y[t, i] * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                                         for t in range(k, t_max[k, i]+1))
                                        >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, t_max[k, i]+1)))
                    if l > 0 and k > 0:
                        left = SR[l-1, j].X/alpha[i][j] + SS[k-1, i].X + sum(Y[t, i].X * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                                                       for t in range(k, t_max[k, i]+1))

                        right = sum((d[i][t] - L[t, i].X)*U[k][t][i]
                                    for t in range(k, t_max[k, i]+1))

                        if left + 0.01 < right:
                            m._ContCut += 1
                            model.addConstr(SR[l-1, j]/alpha[i][j] + SS[k-1, i] + sum(Y[t, i] * min(sum(r[j][p]/alpha[i][j] for p in range(l, t+1)), sum(d[i][v]*U[k][v][i] for v in range(t, t_max[k, i]+1)))
                                                                      for t in range(k, t_max[k, i]+1))
                                        >= sum((d[i][t] - L[t, i])*U[k][t][i] for t in range(k, t_max[k, i]+1)))
                                        
           
def Proposition3(model):
        
            for i in item:
                for j in item:
                    for k in period[0:len(period)]:
                        N = {}
                        left = 0
                        for t in range(k+1, len(period)):
                            if Y[t, i].X * sum(r[j][p] for p in range(k+1, t+1)) < alpha[i][j]*X[t, i].X:
                                N[i, j, t] = 1
                                left += Y[t, i].X * sum(r[j][p] for p in range(k+1, t+1)) - X[t, i].X
                            else:
                                N[i, j, t] = 0

                        if SR[k, j].X + left + 0.01 < 0:
                            m._ContCut += 1
                            model.addConstr(SR[k, j] + sum(N[i, j, t]*Y[t, i] * sum(r[j][p] for p in range(k+1, t+1)) for t in range(k+1, len(period)))
                                        >= sum(alpha[i][j]*X[t, i]*N[i, j, t] for t in range(k+1, len(period))))
            
            
    

# %%INSTANCIAS
TimeHorizon = [25]
product = [5,10]
SetupCostRatio = [1000]
ProductionCostRatio = [2,4]
rdRatio = [1,2,3]
orden = 0

CutGeneration = 0
NCuts=0
RootNode = True
bigMcheck = True
# PARA CALLBACK
da = pd.DataFrame(columns=["instancia", "Time horizon", "Item", "SetupRatio", "ProdRatio", "rdRadio",
                           "MIPGAP", "Tiempo", "TiempoCallback", "ValorObjetivo", "ValorRelax", "Best Bound", "SolucionesOptimas", "Cortes", "Nodos"])

for periodo in TimeHorizon:
    for setupRatio in SetupCostRatio:
        for prodRatio in ProductionCostRatio:
            for rd in rdRatio:
                for itemo in product:
                    z = 1
                    while z <=5:
                        np.random.seed(z*345+z*z)
                        period = range(periodo)
                        item = range(itemo)
                        d = np.round(np.random.uniform(0, 100, (len(item), len(period))),0)
                        # print("d", d)
                        
                        # for ii in item:
                        #     for jj in period:
                        #         if d[ii, jj] < 0:
                        #             d[ii, jj] = 0

                        meanDemand = np.mean(d)
                        r = np.round(np.random.uniform(
                            rd*0.8*meanDemand, rd*1.2*meanDemand, (itemo, periodo)),0)
                        for ii in item:
                            for jj in period:
                                if r[ii, jj] < 0:
                                    r[ii, jj] = 0
                        # print("r", r)
                        # %%PARAMETERS
                        change = np.random.uniform(
                            0.8*setupRatio, 1.2*setupRatio, (len(item), len(item), len(period)))
                        # Initial inventory
                        is_0 = 0
                        ir_0 = 0

                        # Big m
                        mm = 100000

                        p = 200

                        # unit inventory serviceable holding cost in period t.This is charged for inventory at the end of period (h_t)
                        hs_t = np.random.uniform(7, 12, (len(item), len(period)))
                        # unit production (procurement) cost (variable cost) in period t (c_t)
                        ########PARA LA VERSIÓN JOURNAL AGREGAR H BARRITA######
                        c_t = np.random.uniform(0.8*prodRatio, 1.2*prodRatio, (len(item), len(period)))
                        # setup cost (fixed cost) of production(procurement) in period t (s_t)
                        s_t = np.random.uniform(0.8*setupRatio, 1.2*setupRatio, (len(item), len(period)))
                        # Unit inventory remanufacturing holding cost in period t.
                        hr_t = np.random.uniform(2, 7, (len(item), len(period)))
                        # Proportion part of end-of-life product j in the finished product i
                        # proportion part of end-of-life product j in the finished product i
                        alpha = np.random.uniform(1, 6, (itemo, itemo))
                        # for i in item:
                        #     for j in item:
                        #         alpha[i][j] = np.round(1/(len(item)), 2)
                        #         # alpha[i][j] = 1
                        phi = {}
                        for i in item:
                            for k in range(len(period)):
                                # print('\n sum d[%s, %s] : %s'%(k,i,sum(d[i][v]
                                                # for v in range(k, len(period)))))
                                phi[k, i] = sum(d[i][v]
                                                for v in range(k, len(period)))
                                for j in item:
                                    # print('\n sum r[%s, %s, %s] : %s'%(k,i,j,sum(r[j][l]/alpha[i][j] for l in range(k+1))))
                                    
                                    if phi[k, i] > sum(r[j][l]/alpha[i][j] for l in range(k+1)):
                                        phi[k, i] = sum(r[j][l]/alpha[i][j]
                                                        for l in range(k+1))
                        # print("d ", d)                
                        # print("r ", r)
                        # print("alpha ", alpha)
                        # print("phi ",phi)                        
                        # %% VARIABLES
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
                                      vtype=GRB.CONTINUOUS, name="Decision")

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
                                      vtype=GRB.CONTINUOUS, name="Changeover")

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
                        if bigMcheck:
                            m.addConstrs((X[t, i] <= phi[t,i]*Y[t, i]
                                          for t in period for i in item), name="Logical Constraints")
                        else:    
                            m.addConstrs((X[t, i] <= mm*Y[t, i]
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
                        # %%
                        
                        
                        m._time = 0
                        t_modelo = time.time()
                        m._ContCallback = 0
                        m._ContCut = 0
                        m._BestBound = 0
                        m._RootNode = RootNode
                        # VER EL K+1 DE LAS LKU
                        
                        m.setParam("OutputFlag", 0)
                        m.optimize()
                        value = m.objVal
                        checkValue = 0
                        
                        iterations = 0
                        t1 = time.time()
                        while  value - checkValue > 0.1 and iterations <=100:
                            iterations+=1
                            checkValue = value
                            
                            
                            #SingleItem(m)   
                            #SingleItemAlpha(m)
                            #SingleItemReturnsa(m)
                            #Proposition3(m)
                            
                            m.optimize()
                            value = m.objVal
                            print ("#########", iterations, value, m._ContCut)
                        m._time += time.time() - t1 
                        
                            
                        if m.Status == GRB.OPTIMAL:
                            print('\n',i, z, ' Optimal objective: %g' %
                                  m.ObjVal, '\n')

                            tiempo = m.runtime
                            value = m.objVal
                            gap = 0
                            ValueRelax = m.objVal
                            NumSoluciones = m.getAttr("SolCount")
                            Iteraciones = m.getAttr("IterCount")
                            cortes = m.getAttr("NodeCount")
                            #m._BestBound = m.getAttr("ObjBound")
                            # PARA CALLBACK
                            da.loc[orden] = (z, periodo, itemo, setupRatio, prodRatio, rd, gap, tiempo,
                                             m._time, value, ValueRelax, m._BestBound, NumSoluciones, m._ContCut, m._ContCallback)

                        elif m.Status == GRB.INF_OR_UNBD:
                            print('Model is infeasible or unbounded')
                            da.loc[orden] = ("infeasible", "infeasible", "infeasible", 'infeasible', 'infeasible', 'infeasible',
                                             'infeasible', 'infeasible', 'infeasible', 'infeasible', 'infeasible', "q", "q", "q", "q")
                        elif m.Status == GRB.INFEASIBLE:
                            da.loc[orden] = ("infeasible", "infeasible", "infeasible", 'infeasible', 'infeasible', 'infeasible',
                                             'infeasible', 'infeasible', 'infeasible', 'infeasible', 'infeasible', "q", "q", "q", "q")

                        elif m.Status == GRB.UNBOUNDED:
                            da.loc[orden] = ("infeasible", "infeasible", "infeasible", 'infeasible', 'infeasible', 'infeasible',
                                             'infeasible', 'infeasible', 'infeasible', 'infeasible', 'infeasible', "q", "q", "q", "q")

                        else:
                            gap = m.MIPgap
                            tiempo = m.runtime
                            value = m.objVal
                            r = m.relax()
                            r.optimize()
                            ValueRelax = r.ObjVal
                            NumSoluciones = m.getAttr("SolCount")
                            Iteraciones = m.getAttr("IterCount")
                            cortes = m.getAttr("NodeCount")
                            #m._BestBound = m.getAttr("ObjBound")
                            print('Optimization ended with status %d' %
                                  m.Status)
                            da.loc[orden] = (z, periodo, itemo, setupRatio, prodRatio, rd, gap, tiempo,
                                             m._time, value, ValueRelax, m._BestBound, NumSoluciones, m._ContCut, m._ContCallback)

                        orden += 1
                        z += 1

                        da.to_excel("WITHOUTCUT_TIEMPOGUROBI.xlsx")
