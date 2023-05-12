import os, sys
import pyqtgraph as pg
import pyqtgraph.console

from aisynphys.ui.multipatch_nwb_viewer import MultipatchNwbViewer


if __name__ == '__main__':
    # open a console for debugging
    message = """Console variables:

    analyzer : The currently active analyzer
    sweeps : A list of the currently selected sweeps
    view : The MiesNwbViewer instance
    nwb : The currently loaded NWB file

    """
    pg.mkQApp()
    
    console = pg.console.ConsoleWidget(text=message)
    console.catchNextException()
    console.show()

    # create NWB viewer
    v = MultipatchNwbViewer()
    v.show()
    v.setWindowTitle('Multipatch NWB Viewer')

    # set up a convenient function for loading nwb from filename
    nwb = None
    def load_nwb(filename):
        global nwb
        print("Loading %s" % filename)
        nwb = v.load_nwb(filename)
        console.localNamespace['nwb'] = nwb

    # make a few variables available from the console
    console.localNamespace.update({'view': v, 'nwb': None, 'sweeps': [], 'analyzer': None})

    # make selected sweeps / analyzer available as well
    def update_namespace():
        console.localNamespace.update({'sweeps': v.selected_sweeps(), 'analyzer': v.selected_analyzer()})

    v.analyzer_changed.connect(update_namespace)
    v.explorer.selection_changed.connect(update_namespace)

    # load file or set base directory from argv
    for arg in sys.argv[1:]:
        if not os.path.exists(arg):
            from aisynphys.database import default_db as db
            try:
                expt = db.experiment_from_ext_id(arg)
            except Exception:
                print("Could not find experiment %s" % arg)
                sys.exit(-1)
            arg = os.path.join(expt.nwb_file)

        if os.path.isfile(arg):
            load_nwb(arg)

    # start Qt event loop if this is not an interactive python session
    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec_()
