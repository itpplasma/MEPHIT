from os import path
from numpy import pi
from h5py import File
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import magdifplot

workdir = '/home/patrick/git/NEO-EQ/run/33353_2325'
data = File(path.join(workdir, 'magdif.h5'), 'r')
r = data['/cache/fs/rad'][()]
A = data['/cache/fs/area'][()]
C = data['/cache/fs/perimeter'][()]
half_r = data['/cache/fs_half/rad'][()]
half_A = data['/cache/fs_half/area'][()]
half_C = data['/cache/fs_half/perimeter'][()]

fig = Figure()
ax = fig.subplots()
ax.plot(r, A, '-k')
ax.plot(half_r, half_A, '--r')
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$A$ / \si{\centi\meter\squared}')
canvas = FigureCanvas(fig)
fig.savefig(path.join(workdir, 'fs_area.pdf'))

fig = Figure()
ax = fig.subplots()
ax.plot(r, C, '-k')
ax.plot(half_r, half_C, '--r')
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$C$ / \si{\centi\meter}')
canvas = FigureCanvas(fig)
fig.savefig(path.join(workdir, 'fs_perimeter.pdf'))

fig = Figure()
ax = fig.subplots()
ax.plot(r[1:], r[1:] * C[1:] / A[1:], '-k')
ax.plot(half_r, half_r * half_C / half_A, '--r')
ax.hlines([4.0 / pi, 2.0], 0, 1, transform=ax.get_yaxis_transform(), ls=':')
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$r C / A$')
canvas = FigureCanvas(fig)
fig.savefig(path.join(workdir, 'fs_ratio.pdf'))
