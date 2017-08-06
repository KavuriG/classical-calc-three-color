#!usr/bin/python -tt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as col

FUND_FREQ = 3.7474057E14
PL_FWHM = 0
ITER = 0

SOL_LIST = []
t_grid = np.linspace(-5*PL_FWHM, -5*PL_FWHM + 20/FUND_FREQ, 20000)

for t in t_grid[:np.argmax(t_grid >= 1/FUND_FREQ):10]:
    ITER += 1
    a = -1*np.cos(2*np.pi*FUND_FREQ*t)
    b = -1*np.sin(2*np.pi*FUND_FREQ*t)/(2*np.pi*FUND_FREQ) - a*t
    ind = np.argmax(t_grid > t)
    x = np.zeros(ind)
    x[:] = np.NaN
    sol = np.concatenate([x, np.sqrt((np.sin(2*np.pi*FUND_FREQ*t_grid[ind:])/(2*np.pi*FUND_FREQ)
                                      + a*t_grid[ind:] + b)**2)])
    #sol = np.concatenate([x, -1*np.cos(FUND_FREQ*2*np.pi*t_grid[ind:]) + a])
    SOL_LIST.append(sol)



fig1 = plt.figure(1)
ax1 = plt.axes([0.05, 0.15, 0.9, 0.80])
ax1.clear()
cmap = col.copper
for i in range(len(SOL_LIST)):
    l = ax1.plot(t_grid, SOL_LIST[i], color=cmap(i/float(len(SOL_LIST))))
y_field = np.sin(2*np.pi*FUND_FREQ*t_grid)/FUND_FREQ
l1 = ax1.plot(t_grid, y_field, 'b')
#l2 = ax1.plot(time, kin_en, 'g')
#ax1.set_ylim(-0.1E-7, 0.5E-7)
#ax1.set_xlim(-5*PL_FWHM , -5*PL_FWHM + 2/(FUND_FREQ))
print(ITER)
fig1.canvas.draw_idle()
plt.show()
