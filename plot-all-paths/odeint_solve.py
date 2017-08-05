#!/usr/bin/python

import numpy as np
from scipy.integrate import odeint
import sys
import matplotlib.pyplot as plt

ELEC_MASS = 9.10938356E-31
ELEC_CHARGE = -1.60217662E-19
FUND_FREQ = 3.7474057E14
SP_LIGHT = 3E8
FOCUS_RADIUS = 30E-6
PL_FWHM = 5E-15
PULSE_ENERGY = 0.6E-3
EPSILON_o = 8.85418782E-12
TIME_GRID = 20000


INTENSITY = 1.88*(PULSE_ENERGY/(FOCUS_RADIUS**2*PL_FWHM))/np.pi #gaussian
FIELD_AMP = np.sqrt(2*INTENSITY/(EPSILON_o*SP_LIGHT))
FIELD_TOLERANCE = FIELD_AMP*1E-3


def chop_start_end(vector, closeness):
    b = np.where(np.diff(vector < closeness))[0] 
    '''b is now an array of all zero crossings of vector'''
    if len(b) <= 1:
        '''must be a direct electron (vector==distance)'''
        ret = np.zeros(len(vector))
        ret[:] = np.Inf
        return ret
    else:
        b_start = b[1]
        '''first zero crossing always AT 0, so choose second crossing'''
        diff = np.diff(vector[b_start:])
        a = np.where(np.diff(np.sign(diff)))[0]
        '''a contains locations of turning points after b_start'''
        if len(a) <= 1:
            ret = np.zeros(len(vector))
            ret[:] = np.Inf
            return ret
        else:
            a_end = a[1]
            '''first turning point is the closest distance,
            choose second to have a good looking plot'''
#            print(a_end,b_start, diff)
            front = vector[:b_start]
            back = vector[b_start + a_end:]
            front[:] = np.NaN
            back[:] = np.NaN
            return np.concatenate([front, vector[b_start:b_start + a_end], back])
 

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
            hard_end = start_i + 20/(FUND_FREQ)
        interval = (hard_end - start_i)/(TIME_GRID - 1)
        t = np.linspace(start_i, hard_end, TIME_GRID)
        t_chopped = t[:np.argmax(t >= end_i)]
        #print(start_i, end_i, hard_end, (t > end_i))
        sol_y = odeint(elec_path_y, initial_cond, t)
        sol_z = odeint(elec_path_z, initial_cond, t)
        plt.plot(sol_y[:,1], 'b.')
        throw = np.zeros(TIME_GRID - len(t_chopped))
        throw[:] = np.NaN
        '''these solve the diff eqn'''
        start_e = start_i
        for i in range(len(t_chopped)):
            '''Optimize over multiple start times:
            as the ode soln. is simply integrating the field twice,
            we can easily generate solns. for diff. start times
            by applying diff. boundary conditions and changing
            the constants of integration. B.C. is always 0 pos & 0 vel.
            at the start time. This is better than repeatedly calling odeint'''
            a_y = -1*sol_y[i, 1]
            if i==0:
                print(sol_y, i)
            b_y = -1*(sol_y[i, 0] + a_y*t[i])
            a_z = -1*sol_z[i, 1]
            b_z = -1*(sol_z[i, 0] + a_z*t[i])
            y = sol_y[:, 0] + a_y*t + b_y
            z = sol_z[:, 0] + a_z*t + b_z
            #sol_y[:,1] = sol_y[:,1] + a_y
            #sol_z[:,1] = sol_z[:,1] + a_z
            #end
            if i%10 == 0:
                a = np.zeros(i)
                a[:] = np.NaN
                dist = np.concatenate([a, np.sqrt((y[i:])**2)])
                                                  #+ z[i:]**2)])
                min_list.append((dist, start_e))
            #min_list.append((kin_energy, start_i))
            start_e += interval
        minimas = np.array(min_list)
#        min_start = minimas[np.nanargmin(minimas, axis=0)[0], 1]
        print(start_i, end_i)
        #t = np.linspace(start_i, end_i, TIME_GRID)
#        sol_y = odeint(elec_path_y, initial_cond, t, full_output=0)
#        sol_z = odeint(elec_path_z, initial_cond, t, full_output=0)
        return [minimas, t]
#        return minimas
    else:
        t = np.linspace(start_i, hard_end, TIME_GRID)
        sol_y = odeint(elec_path_y, initial_cond, t)
        sol_z = odeint(elec_path_z, initial_cond, t)
return [0, sol_y[:,0], sol_z[:,0]]
