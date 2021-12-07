import h5py
from os import path
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from magdifplot import run_dir
import numpy as np

work_dir = run_dir + '/33353_2325'
data = h5py.File(path.join(work_dir, 'mephit.h5'), 'r')
m_res_min = data['/mesh/m_res_min'][()]
m_res_max = data['/mesh/m_res_max'][()]
res_ind = data['/mesh/res_ind'][()] - 1
kt_max = data['/mesh/kt_max'][()]
kt_low = data['/mesh/kt_low'][()] - 1
npoint = data['/mesh/npoint'][()]
GL_order = data['/mesh/GL_order'][()]
GL_R = data['/mesh/GL_R'][()]
GL_Z = data['/mesh/GL_Z'][()]
shielding_kt_low = data['/mesh/shielding_kt_low'][()] - 1
shielding_L1_R = data['/mesh/shielding_L1_R'][()]
shielding_L1_Z = data['/mesh/shielding_L1_Z'][()]
fig = Figure(figsize=(3.3, 4.4))
ax = fig.subplots()
ax.triplot(data['/mesh/node_R'][()], data['/mesh/node_Z'][()], data['/mesh/tri_node'][()] - 1, lw=0.01)
for m in range(m_res_min, m_res_max + 1):
    km = m - m_res_min
    kf = res_ind[km]
    for kt in range(kt_max[kf]):
        kedge = npoint + kt_low[kf] + kt
        k = shielding_kt_low[km] + kt + 1
        for kgl in range(GL_order):
            R = np.array([shielding_L1_R[k, kgl, 0], GL_R[kedge, kgl], shielding_L1_R[k, kgl, 1]])
            Z = np.array([shielding_L1_Z[k, kgl, 0], GL_Z[kedge, kgl], shielding_L1_Z[k, kgl, 1]])
            ax.plot(R, Z, '-r', lw=0.01)
ax.set_aspect('equal')
ax.set_xlabel(r'$R$ / cm')
ax.set_ylabel(r'$Z$ / cm')
canvas = FigureCanvas(fig)
fig.savefig(work_dir + '/shielding.pdf', dpi=300)
