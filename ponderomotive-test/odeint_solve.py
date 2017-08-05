#!/usr/bin/python

import numpy as np
from scipy.integrate import odeint
import sys

ELEC_MASS = 9.10938356E-31
ELEC_CHARGE = -1.60217662E-19
FUND_FREQ = 3.7474057E14
SP_LIGHT = 3E8
FOCUS_RADIUS = 30E-6
PL_FWHM = 5E-15
PULSE_ENERGY = 0.6E-3
EPSILON_o = 8.85418782E-12
TIME_GRID = 5000


INTENSITY = 1.88*(PULSE_ENERGY/(FOCUS_RADIUS**2*PL_FWHM))/np.pi #gaussian
FIELD_AMP = np.sqrt(2*INTENSITY/(EPSILON_o*SP_LIGHT))
FIELD_TOLERANCE = FIELD_AMP*1E-3


def chop_start_end(vector, closeness):
    '''cut the given vector to the foll. specs:
        FROM the second time the value is less than closeness
        TO the second point of zero derivative after that.
        
        Used to chop trajectories at their first collision with the
        atom, but leave large room at the ends so that it can be 
        plotted if need be.'''
    b = np.where(np.diff(vector < closeness))[0]
    if len(b) <= 1:
        return np.array([np.inf])
    else:
        b_start = b[1]
        diff = np.ediff1d(vector[b_start:])
        a = np.where(np.diff(np.sign(diff)))[0]
        if len(a) <= 1:
            return np.array([np.inf])
        else:
            a_end = a[1]
#            print(a_end,b_start, diff)
            return vector[b_start:b_start+a_end]


def find_end(field, start, end, tol):
    interval = (end - start)*10/TIME_GRID
    probe = end
    while np.sqrt((field[0](probe))**2 + (field[1](probe))**2) > tol:
        probe+= interval
    return probe

def solve_path(field, start_i, end_i, optimize_collision=True,
               pulsed=True, closeness=1, *args):
    '''Solves the eqn for the electron path, returns closest approach

    Args:
        field (list, [y_comp, z_comp]),
        start_ionize_time, end_ionize time,
        bool (optimize collision?), bool (pulsed?)
    Args(optional):
        args[0] = cycle length of CW light
        (compulsory if optimize_collision and not pulsed)
        TODO: change time grid as a function input'''

    def elec_path_y(x, t):
           path, velocity = x
           dydt = [velocity, ELEC_CHARGE*field[0](t)/ELEC_MASS]
           return dydt
    def elec_path_z(x, t):
           path, velocity = x
           dzdt = [velocity, ELEC_CHARGE*field[1](t)/ELEC_MASS]
           return dzdt
    initial_cond = [0., 0.]
    if optimize_collision:
        hard_end = 0
        min_list = []
        if pulsed:
            hard_end = find_end(field, start_i, end_i, FIELD_TOLERANCE)
        else:
            #hard_end = 6*PL_FWHM
            hard_end = start_i + 20.0/(FUND_FREQ)
        interval = (hard_end - start_i)/(TIME_GRID - 1)
        t = np.linspace(start_i, hard_end, TIME_GRID)
        t_chopped = t[:np.argmax(t >= end_i)]
        #print(start_i, end_i, hard_end, (t > end_i))
        sol_y = odeint(elec_path_y, initial_cond, t)
        sol_z = odeint(elec_path_z, initial_cond, t)
        start_e = start_i        
        for i in range(len(t_chopped)):
            '''Optimize over multiple start times'''
            a_y = -1*sol_y[i,1]
            b_y = -1*(sol_y[i,0] + a_y*t[i])
            a_z = -1*sol_z[i,1]
            b_z = -1*(sol_z[i,0] + a_z*t[i])
            y_pos = sol_y[:,0] + a_y*t + b_y
            z_pos = sol_z[:,0] + a_z*t + b_z
            y_vel = sol_y[:,1] + a_y
            z_vel = sol_z[:,1] + a_z            
            #end
            #dist = np.sqrt(sol_y[i:, 0]**2 + sol_z[i:, 0]**2)
            if len(np.where(np.diff(np.sign(y_pos[i:])))[0]) > -1:
                #print("Whaaa!!", sol_y[i:,0])
                #ret_ind = np.where(np.diff(np.sign(sol_y[i:,0])))[0][1]
                kin_energy = 0.5*ELEC_MASS*(y_vel[-1]**2 + z_vel[-1]**2)
                #print(kin_energy)
                '''dist = chop_start_end(np.sqrt(sol_y[i:, 0]**2 + sol_z[i:, 0]**2), closeness)'''
                '''min_list.append((np.amin(dist), start_i, kin_energy))'''
                min_list.append((kin_energy, start_e))
            start_e += interval
        minimas = np.array(min_list)
#        min_start = minimas[np.argmin(minimas, axis=0)[0], 1]
#        t = np.linspace(min_start, hard_end, TIME_GRID)            
#        sol_y = odeint(elec_path_y, initial_cond, t, full_output=0)
#        sol_z = odeint(elec_path_z, initial_cond, t, full_output=0)
#        dist = chop_start_end(np.sqrt(sol_y[:, 0]**2 + sol_z[:, 0]**2), closeness)
#        dis2 = np.sqrt(sol_y[:, 0]**2 + sol_z[:, 0]**2)
#        return [min_start, sol_y[:,0], sol_z[:,0],dist, dis2]
        return minimas
    else:
        t = np.linspace(start_i, end_i, TIME_GRID)            
        sol_y = odeint(elec_path_y, initial_cond, t)
        sol_z = odeint(elec_path_z, initial_cond, t)
        return [0, sol_y[:,0], sol_z[:,0]]