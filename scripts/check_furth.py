from magdifplot import magdif
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from os import path
import numpy as np

work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'
testcase = magdif(work_dir)
testcase.read_datafile()
rad = testcase.data['/debug_furth/rad'][()]
k_z = testcase.data['/debug_furth/k_z'][()]
k_theta = testcase.data['/debug_furth/k_theta'][()]
Bmn_plas_rad = testcase.data['/debug_furth/Bmn_plas_rad'][()]
furth_2 = testcase.data['/debug_furth/sheet_flux'][()]
I = np.delete(furth_2, furth_2 == 0.0 + 0.0j)
k2 = k_theta * k_theta + k_z * k_z
furth_1 = np.empty(shape=furth_2.shape, dtype=complex)
furth_1.real = rad / k2 * np.gradient(Bmn_plas_rad.real, rad)
furth_1.imag = rad / k2 * np.gradient(Bmn_plas_rad.imag, rad)
kmax_Re = np.argmax(furth_1.real)
kmax_Im = np.argmax(furth_1.imag)
kmin_Re = np.argmin(furth_1.real)
kmin_Im = np.argmin(furth_1.imag)
jump_Re = (furth_1[kmax_Re].real - furth_1[kmin_Re].real) * np.sign(kmax_Re - kmin_Re)
jump_Im = (furth_1[kmax_Im].imag - furth_1[kmin_Im].imag) * np.sign(kmax_Im - kmin_Im)
print(f"jump_Re = {jump_Re}, jump_Im = {jump_Im}")
fig = Figure()
ax = fig.subplots()
ax.plot(rad, furth_1.real, '-k', label='real')
ax.plot(rad, furth_1.imag, '--r', label='imag')
ax.plot(rad[[kmax_Re, kmin_Re]], furth_1[[kmax_Re, kmin_Re]].real, 'xk')
ax.plot(rad[[kmax_Im, kmin_Im]], furth_1[[kmax_Im, kmin_Im]].imag, 'xr')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, "check_furth.pdf"))
