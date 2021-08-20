from os import path
from numpy import pi, sqrt, zeros
from h5py import File
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import magdifplot

workdir = '/home/patrick/git/NEO-EQ/run/33353_2325'
data = File(path.join(workdir, 'magdif_fix.h5'), 'r')
nflux = data['/mesh/nflux'][()]
rad = data['/cache/fs/rad'][()]
A = data['/cache/fs/area'][()]
C = data['/cache/fs/perimeter'][()]
half_rad = data['/cache/fs_half/rad'][()]
half_A = data['/cache/fs_half/area'][()]
half_C = data['/cache/fs_half/perimeter'][()]

Delta_r = zeros((nflux))
for kf in range(1, nflux):
    Delta_r[kf] = (A[kf] - A[kf - 1]) / half_C[kf - 1]
Delta_half_r = zeros((nflux))
Delta_half_r[0] = 2.0 * half_A[0] / half_C[0]
for kf in range(1, nflux - 1):
    Delta_half_r[kf] = (half_A[kf] - half_A[kf - 1]) / C[kf]
Delta_half_r[nflux - 1] = 2.0 * (A[nflux] - half_A[nflux - 1]) / C[nflux]
var_half_r = sqrt(half_A / pi)
Delta_var_half_r = var_half_r[1:] - var_half_r[:-1]

fig = Figure()
ax = fig.subplots()
ax.plot(rad, A, '-k')
ax.plot(half_rad, half_A, '--r')
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$A$ / \si{\centi\meter\squared}')
canvas = FigureCanvas(fig)
fig.savefig(path.join(workdir, 'fs_area.pdf'))

fig = Figure()
ax = fig.subplots()
ax.plot(rad, C, '-k')
ax.plot(half_rad, half_C, '--r')
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$C$ / \si{\centi\meter}')
canvas = FigureCanvas(fig)
fig.savefig(path.join(workdir, 'fs_perimeter.pdf'))

fig = Figure()
ax = fig.subplots()
ax.plot(rad[1:], rad[1:] * C[1:] / A[1:], '-k')
ax.plot(half_rad, half_rad * half_C / half_A, '--r')
ax.hlines([4.0 / pi, 2.0], 0, 1, transform=ax.get_yaxis_transform(), ls=':')
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$r C / A$')
canvas = FigureCanvas(fig)
fig.savefig(path.join(workdir, 'fs_ratio.pdf'))

fig = Figure()
ax = fig.subplots()
ax.plot(half_rad[1:], rad[2:] - rad[1:-1], '-k')
ax.plot(rad[1:], Delta_half_r * 2.0 / sqrt(3.0), '--r')
ax.plot(rad[1:-1], (Delta_r[1:] + Delta_r[:-1]) / sqrt(3.0), '-.b')
ax.plot(rad[1:-1], Delta_var_half_r * 2.0 / sqrt(3.0), ':g')
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$l_{k}$ / \si{\centi\meter}')
canvas = FigureCanvas(fig)
fig.savefig(path.join(workdir, 'optpolres.pdf'))
