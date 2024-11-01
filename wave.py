import math
import matplotlib.pyplot as plt
import numpy as np

N = 2**7
L = 2 * math.pi
h = L / N
x = [h*e for e in range(1, N + 1)]
c = [0.2 + math.sin(e-1)**2 for e in x]
dt = h / 4

tmax = 16
tplot = 0.15
plotgap = round(tplot/dt)
dt = tplot/plotgap
nplots = round(tmax/tplot)

v = [math.exp(-10*(e-1)**2) for e in x]
vold = [math.exp(-10*(e-0.2*dt-1)**2) for e in x]

k = range(N//2 + 1)
k1 = [1j*e for e in k]; k1[-1] = 0
k2 = [-e**2 for e in k]

t = 0
l, = plt.plot(x, v, 'o-')
plt.xlim([0, 2*math.pi])
plt.ylim([-0.2, 1.2])
for i in range(nplots):
    for n in range(plotgap):
        t += dt
        v_hat = np.fft.rfft(v)
        w_hat = v_hat * k1
        w = np.fft.irfft(w_hat)
        vnew = vold - 2*dt*w*c
        vold = v
        v = vnew
    l.set_xdata(x)
    l.set_ydata(v)
    plt.savefig("u.%08d.png" % i)
