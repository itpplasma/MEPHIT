from mephit_plot import Gpec, Mephit, ParallelPlotter, Plot1D, \
    run_dir, set_matplotlib_defaults, XTicks
from numpy import sign


if __name__ == "__main__":
    set_matplotlib_defaults()
    plotter = ParallelPlotter()
    plotter.start()

    work_dir = run_dir + '/33353_2900_EQH'
    testcase = Mephit(work_dir)
    testcase.open_datafile()
    testcase.postprocess()
    sgn_dpsi = sign(testcase.data['/cache/fs/psi'][-1] - testcase.data['/cache/fs/psi'][0])
    nflux = testcase.data['/mesh/nflux'][()]

    mephit_Ires = testcase.get_Ires()
    reference = Gpec(work_dir, 2)
    reference.open_datafiles()
    gpec_Ires = reference.get_Ires()
    config = {
        'xlabel': '$m$', 'ylabel': r'$\abs\, I_{m, n}^{\parallel}$ / \si{\ampere}', 'legend': {'fontsize': 'small'},
        'plotdata': [
            {'x': mephit_Ires.keys(), 'y': mephit_Ires.values(), 'args': {'label': 'MEPHIT', 'marker': 'o', 'ls': ''}},
            {'x': gpec_Ires.keys(), 'y': gpec_Ires.values(), 'args': {'label': 'GPEC', 'marker': 'x', 'ls': ''}}
        ],
        'postprocess': [XTicks(mephit_Ires.keys())]
    }
    plotter.plot_objects.put(Plot1D(work_dir, 'plot_Ires.pdf', config))

    plotter.finish()
    testcase.close_datafile()
    reference.close_datafiles()
