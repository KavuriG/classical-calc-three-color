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

PONDER = (FIELD_AMP**2)*(ELEC_CHARGE**2)/(4*ELEC_MASS*(2*np.pi*FUND_FREQ)**2)


def plot(*args, **kwargs):
#    x = np.linspace(slider_12.val*PL_FWHM,(slider_12.val + slider_11.val)*PL_FWHM,200)

#    z_field = args[1][1](x)
    t2 = np.linspace(args[1], args[2], len(args[-1]))
#    y = args[0][0]
#    z = args[0][1]
    dist = args[-1]
#    closest = 1E-9*args[0][:,0]/PONDER
#    mask = [0 if np.ma.is_masked(i) else np.NaN for i in np.ma.masked_invalid(closest)]
#    time = args[0][:,1]
    #kin_en = 1E-9*args[0][:,2]/np.max(args[0][:,2])
    y_field = 1E-9*kwargs['field'][0](t2)/FIELD_AMP
    ax1.clear()
#    l = ax1.plot(time, closest, 'b')
#    l2 = ax1.plot(time, mask, 'r.')
    l = ax1.plot(t2, dist, 'b',t2, args[-2], 'r')
    l1 = ax1.plot(t2,y_field, 'g')
    #l2 = ax1.plot(time, kin_en, 'g')
    ax1.set_ylim(-1E-9, 2E-9)
    ax1.set_xlim(-5*PL_FWHM, -4*PL_FWHM)
    fig1.canvas.draw_idle()




def update(val):

    qwp_1 = slider_1.val
    hwp_2 = slider_2.val
    qwp_2 = slider_3.val
    hwp_3 = slider_4.val
    qwp_3 = slider_5.val
    delay_1 = slider_6.val
    delay_2 = slider_7.val
    ampl_1 = slider_8.val
    ampl_2 = slider_9.val
    ampl_3 = slider_10.val
    closeness = slider_13.val

    a = e_field.e_field_gen(3, False, ampl_1*FIELD_AMP, ampl_2*FIELD_AMP,
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
    times = t#[j for i,j in zip(tot,t) if i > FIELD_AMP_ION]
    b = ode.solve_path(a, times[0], times[0] + 2/(FUND_FREQ), True, False, closeness*1E-9)
    #print(b)
#    plot(b,a)
    plot(b[2:4], b[0], b[1], b[4], b[5], field=a)
#    plot(a)


if __name__ == '__main__':


    fig1 = plt.figure(1)
    ax1 = plt.axes([0.05, 0.15, 0.9, 0.80])
    #fig2 = plt.figure(2)
    #ax2 = plt.axes([0.05, 0.15, 0.9, 0.80], projection='3d')
    ax_slider_1 = plt.axes([0.1, 0.01, 0.2, 0.02])
    ax_slider_2 = plt.axes([0.1, 0.04, 0.2, 0.02])
    ax_slider_3 = plt.axes([0.1, 0.07, 0.2, 0.02])
    ax_slider_4 = plt.axes([0.1, 0.1, 0.2, 0.02])
    ax_slider_5 = plt.axes([0.1, 0.13, 0.2, 0.02])
    ax_slider_6 = plt.axes([0.5, 0.01, 0.2, 0.02])
    ax_slider_7 = plt.axes([0.5, 0.04, 0.2, 0.02])
    ax_slider_8 = plt.axes([0.5, 0.07, 0.2, 0.02])
    ax_slider_9 = plt.axes([0.5, 0.1, 0.2, 0.02])
    ax_slider_10 = plt.axes([0.5, 0.13, 0.2, 0.02])
    rax = plt.axes([0.0, 0.8, 0.1, 0.15])
    ax_slider_11 = plt.axes([0.2, 0.97, 0.1, 0.02])
    ax_slider_12 = plt.axes([0.6, 0.97, 0.1, 0.02])
    ax_slider_13 = plt.axes([0.8, 0.97, 0.2, 0.02])
    radio = widgets.RadioButtons(rax, ('CW', 'Pulsed'))
    slider_1 = widgets.Slider(ax_slider_1, 'qwp_1', 0., 360)
    slider_2 = widgets.Slider(ax_slider_2, 'hwp_2', 0., 360)
    slider_3 = widgets.Slider(ax_slider_3, 'qwp_2', 0., 360)
    slider_4 = widgets.Slider(ax_slider_4, 'hwp_3', 0., 360)
    slider_5 = widgets.Slider(ax_slider_5, 'qwp_3', 0., 360)
    slider_6 = widgets.Slider(ax_slider_6, 'delay_2', -2, 2)
    slider_7 = widgets.Slider(ax_slider_7, 'delay_3', -2, 2)
    slider_8 = widgets.Slider(ax_slider_8, 'ampl_1', 0, 1)
    slider_9 = widgets.Slider(ax_slider_9, 'ampl_2', 0, 1)
    slider_10 = widgets.Slider(ax_slider_10, 'ampl_3', 0, 1)
    slider_11 = widgets.Slider(ax_slider_11, 'x-size', 0, 4)
    slider_12 = widgets.Slider(ax_slider_12, 'x-start', -2, 2)
    slider_13 = widgets.Slider(ax_slider_13, 'close', 0, 1)
    #start


    qwp_1 = slider_1.val
    hwp_2 = slider_2.val
    qwp_2 = slider_3.val
    hwp_3 = slider_4.val
    qwp_3 = slider_5.val
    delay_1 = slider_6.val
    delay_2 = slider_7.val
    ampl_1 = slider_8.val
    ampl_2 = slider_9.val
    ampl_3 = slider_10.val
    closeness = slider_13.val

    a = e_field.e_field_gen(3, False, ampl_1*FIELD_AMP, 0*FIELD_AMP,
                            0*FIELD_AMP, 0,
                            delay_1/FUND_FREQ, delay_2/FUND_FREQ,
                            FUND_FREQ, 2*FUND_FREQ, 3*FUND_FREQ,
                            [[45], [0, 0], [0, 0]],
                            PL_FWHM, PL_FWHM, PL_FWHM,
                            b1='q', b2='hq', b3='hq')

    t = np.linspace(-5*PL_FWHM, 5*PL_FWHM, 50)
    y_field = a[0](t)
    z_field = a[1](t)
    tot = np.sqrt(y_field**2 + z_field**2)
    times = t#[j for i,j in zip(tot,t) if i > FIELD_AMP_ION]
    b = ode.solve_path(a, times[0], times[0] + 2/(FUND_FREQ), True, False, closeness*1E-9)
    #print(b)
#    plot(b,a)
    plot(b[2:4], b[0], b[1], b[4], b[5], field=a)
#    plot(a)



    #end
    slider_1.on_changed(update)
    slider_2.on_changed(update)
    slider_3.on_changed(update)
    slider_4.on_changed(update)
    slider_5.on_changed(update)
    slider_6.on_changed(update)
    slider_7.on_changed(update)
    slider_8.on_changed(update)
    slider_9.on_changed(update)
    slider_10.on_changed(update)
    slider_11.on_changed(update)
    slider_12.on_changed(update)
    slider_13.on_changed(update)

    radio.on_clicked(update)


    plt.show()