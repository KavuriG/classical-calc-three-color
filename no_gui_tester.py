#!usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import e_field_gen as e_field
import odeint_solve as ode
import sys as sys

ELEC_MASS = 9.10938356E-31
ELEC_CHARGE = -1.60217662E-19
FUND_FREQ = 3.7474057E14
SP_LIGHT = 3E8
PL_FWHM = 25E-15
FOCUS_RADIUS = 30E-6
PULSE_ENERGY = 0.6E-3
EPSILON_o = 8.85418782E-12
TIME_GRID = 200

INTENSITY = 1.88*(PULSE_ENERGY/(FOCUS_RADIUS**2*PL_FWHM))/np.pi #gaussian
FIELD_AMP = np.sqrt(2*INTENSITY/(EPSILON_o*SP_LIGHT))
FIELD_TOLERANCE = FIELD_AMP*1E-1
FIELD_AMP_ION = np.sqrt(2E14/(EPSILON_o*SP_LIGHT))


def plot(*args):
#    x = np.linspace(slider_12.val*PL_FWHM,(slider_12.val + slider_11.val)*PL_FWHM,200)
#    y_field = args[0][0](x)
#    z_field = args[0][1](x)
#    t2 = np.linspace(args[1], args[2][-1], len(args[-1]))
#    t = np.linspace(args[1] + 
#                    (len(args[-1]) - len(args[-2]))*(args[2][-1] - args[1])/len(args[-1]),
#                    args[2][-1], len(args[-2]))
#    y = args[0][0]
#    z = args[0][1]
#    dist = args[-1]
        
    fig1 = plt.figure(1)
    ax1 = plt.axes([0.05, 0.15, 0.9, 0.80])
    #fig2 = plt.figure(2)
    #ax2 = plt.axes([0.05, 0.15, 0.9, 0.80], projection='3d')
    
    closest = args[0][:,0]
    mask = [0 if np.ma.is_masked(i) else np.NaN for i in np.ma.masked_invalid(closest)]
    time = args[0][:,1]
    ax1.clear()
    l = ax1.plot(time, closest, 'b.')
    l2 = ax1.plot(time, mask, 'r.')
#    l = ax1.plot(t2, dist, 'r-', t, args[-2], 'b.')
#    l = ax1.plot(x,y_field,z_field)
    ax1.set_ylim(-1E-9, 2E-9)
    ax1.set_xlim(-1E-16, 3E-14)
    plt.savefig(str(np.random.rand(1)[0])+'foo.png')




def update():


    qwp_1 = np.random.rand(1)[0]*360
    hwp_2 = np.random.rand(1)[0]*360
    qwp_2 = np.random.rand(1)[0]*360
    hwp_3 = np.random.rand(1)[0]*360
    qwp_3 = np.random.rand(1)[0]*360
    delay_1 = (np.random.rand(1)[0]-0.5)*4
    delay_2 = (np.random.rand(1)[0]-0.5)*4
    ampl_1 = np.random.rand(1)
    ampl_2 = np.random.rand(1)
    ampl_3 = np.random.rand(1)
    closeness = np.random.rand(1)


    a = e_field.e_field_gen(3, True, ampl_1*FIELD_AMP, ampl_2*FIELD_AMP,
                            ampl_3*FIELD_AMP, 0,
                            delay_1/FUND_FREQ, delay_2/FUND_FREQ,
                            FUND_FREQ, 2*FUND_FREQ, 3*FUND_FREQ,
                            [[qwp_1], [hwp_2, qwp_2], [hwp_3, qwp_3]],
                            PL_FWHM, PL_FWHM, PL_FWHM,
                            b1='q', b2='hq', b3='hq')

    t = np.linspace(-5*PL_FWHM, 5*PL_FWHM, 50)
    y_field = a[0](t)
    z_field = a[1](t)
    tot = np.sqrt(y_field**2 + z_field**2)
    times = [j for i,j in zip(tot,t) if i > FIELD_AMP_ION]
    b = ode.solve_path(a, times[0], times[-1], True, True, closeness*1E-9)
    plot(b)
#    plot(b[1:-2], b[0], times, b[-1], b[-2])
#    plot(a)


if __name__ == '__main__':


    #start


    qwp_1 = np.random.rand(1)[0]*360
    hwp_2 = np.random.rand(1)[0]*360
    qwp_2 = np.random.rand(1)[0]*360
    hwp_3 = np.random.rand(1)[0]*360
    qwp_3 = np.random.rand(1)[0]*360
    delay_1 = (np.random.rand(1)[0]-0.5)*4
    delay_2 = (np.random.rand(1)[0]-0.5)*4
    ampl_1 = np.random.rand(1)
    ampl_2 = np.random.rand(1)
    ampl_3 = np.random.rand(1)
    closeness = np.random.rand(1)

    a = e_field.e_field_gen(3, True, 0*FIELD_AMP, 0*FIELD_AMP,
                            ampl_3*FIELD_AMP, 0,
                            delay_1/FUND_FREQ, delay_2/FUND_FREQ,
                            FUND_FREQ, 2*FUND_FREQ, 3*FUND_FREQ,
                            [[qwp_1], [hwp_2, qwp_2], [hwp_3, qwp_3]],
                            PL_FWHM, PL_FWHM, PL_FWHM,
                            b1='q', b2='hq', b3='hq')

    t = np.linspace(-5*PL_FWHM, 5*PL_FWHM, 50)
    y_field = a[0](t)
    z_field = a[1](t)
    tot = np.sqrt(y_field**2 + z_field**2)
    times = [j for i,j in zip(tot,t) if i > FIELD_AMP_ION]
    b = ode.solve_path(a, times[0], times[-1], True, True, closeness*1E-9)
    plot(b)
#    plot(b[1:-2], b[0], times, b[-1], b[-2])
#    plot(a)



    plt.show()
    
    for i in range(10):
        update()
        
