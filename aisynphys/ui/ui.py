from collections import OrderedDict
import os, sys, subprocess, datetime
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.miesnwb import MiesNwb
from .. import constants
from .. import config


class SynapseTreeWidget(QtWidgets.QTreeWidget):
    """Displays a list of known synapses, given a lit of experiments.
    
    Provides ui for filtering and selecting.
    """
    def __init__(self, expts, parent=None):
        QtWidgets.QTreeWidget.__init__(self, parent)
        self.setColumnCount(3)
        self.expts = expts
        for syn in expts.connection_summary():
            item = SynapseTreeItem(syn['expt'], syn['cells'])
            self.addTopLevelItem(item)

    
class SynapseTreeItem(QtWidgets.QTreeWidgetItem):
    def __init__(self, expt, cells):
        self.expt = expt
        self.cells = cells

        fields = [
            "%0.3f" % expt.timestamp,
            "%d - %d" % (cells[0].cell_id, cells[1].cell_id),
            "%s - %s" % (cells[0].cre_type, cells[1].cre_type),
            
        ]
        
        QtWidgets.QTreeWidgetItem.__init__(self, fields)


class ExperimentInfoWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        
        self.expt = None
        
        self.layout = QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        
        self.info_tree = pg.DataTreeWidget()
        self.layout.addWidget(self.info_tree, self.layout.rowCount(), 0)
        
        self.biocytin_btn = QtWidgets.QPushButton('biocytin image...')
        self.biocytin_btn.clicked.connect(self.show_biocytin)
        self.layout.addWidget(self.biocytin_btn, self.layout.rowCount(), 0)
        
    def set_experiment(self, expt):
        self.expt = expt
        
        info = OrderedDict([
            ('date', str(expt.date)),
            ('specimen', expt.specimen_name),
            ('age', expt.age),
        ])
          
        self.info_tree.setData(info)
        
    def show_biocytin(self):
        if sys.platform == 'win32':
            subprocess.Popen([r'C:\Program Files (x86)\Mozilla Firefox\firefox.exe', self.expt.biocytin_image_url])
        else:
            os.system("firefox " + self.expt.biocytin_image_url)


class ExperimentTimeline(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.channels = None
        self.start_time = None  # starting time according to NWB file
        
        self.layout = QtWidgets.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.plots = PlotGrid()
        self.plots.set_shape(config.n_headstages, 1)
        self.plots.setXLink(self.plots[0, 0])
        self.layout.addWidget(self.plots, 0, 0)

        self.ptree = pg.parametertree.ParameterTree(showHeader=False)
        self.layout.addWidget(self.ptree, 0, 1)
        self.ptree.setMaximumWidth(250)
        
        self.params = pg.parametertree.Parameter.create(name='params',
            type='group', addText='Add pipette')
        self.params.addNew = self.add_pipette_clicked  # monkey!
        self.ptree.setParameters(self.params, showTop=False)
        
    def list_channels(self):
        return self.channels
        
    def get_channel_plot(self, chan):
        return self.plots[chan, 0]
        
    def add_pipette_clicked(self):
        self.add_pipette(channel=self.channels[0], start=0, stop=500)
        
    def remove_pipettes(self):
        for ch in self.params.children():
            self.params.removeChild(ch)
            ch.region.scene().removeItem(ch.region)
        for i in range(self.plots.shape[0]):
            self.plots[i,0].clear()

    def load_site(self, site_dh):
        """Generate pipette list for this site
        """
        self.remove_pipettes()
        
        # automatically fill pipette fluorophore field
        expt_dh = site_dh.parent().parent()
        expt_info = expt_dh.info()
        dye = expt_info.get('internal_dye', None)
        internal = expt_info.get('internal', None)

        # automatically select electrode regions
        self.channels = list(range(config.n_headstages))
        site_info = site_dh.info()
        for i in self.channels:
            headstages = site_info.get('headstages', None)
            if headstages is None:
                hs_state = site_info.get('Headstage %d'%(i+1), None)
            else:
                hs = headstages.get('HS%d'%(i+1), {})
                hs_state = hs.get('Seal', None)
            status = {
                'NS': 'No seal',
                'LS': 'Low seal',
                'GS': 'GOhm seal',
                'TF': 'Technical failure',
                'NA': 'No attempt',
                None: 'Not recorded',                
            }.get(hs_state, hs_state)
            self.add_pipette(i, status=status, internal_dye=dye, internal=internal)
        
    def load_nwb(self, nwb_handle):
        with pg.BusyCursor():
            self._load_nwb(nwb_handle)
        
    def _load_nwb(self, nwb_handle):
        self.nwb_handle = nwb_handle
        self.nwb = MiesNwb(nwb_handle.name())
        
        # load all recordings
        recs = {}
        for srec in self.nwb.contents:
            for chan in srec.devices:
                recs.setdefault(chan, []).append(srec[chan])

        chans = sorted(recs.keys())
        
        # find time of first recording
        start_time = min([rec[0].start_time for rec in recs.values()])
        self.start_time = start_time
        end_time = max([rec[-1].start_time for rec in recs.values()])
        self.plots.setXRange(0, (end_time-start_time).seconds)
        
        # plot all recordings
        for i,chan in enumerate(chans):
            n_recs = len(recs[chan])
            times = np.empty(n_recs)
            i_hold = np.empty(n_recs)
            v_hold = np.empty(n_recs)
            v_noise = np.empty(n_recs)
            i_noise = np.empty(n_recs)
            
            # load QC metrics for all recordings
            for j,rec in enumerate(recs[chan]):
                dt = (rec.start_time - start_time).seconds
                times[j] = dt
                v_hold[j] = rec.baseline_potential
                i_hold[j] = rec.baseline_current
                if rec.clamp_mode == 'vc':
                    v_noise[j] = np.nan
                    i_noise[j] = rec.baseline_noise_stdev
                else:
                    v_noise[j] = rec.baseline_noise_stdev
                    i_noise[j] = np.nan
                    
            # scale all qc metrics to the range 0-1
            pass_brush = pg.mkBrush(100, 100, 255, 200)
            fail_brush = pg.mkBrush(255, 0, 0, 200)
            v_hold = (v_hold + 60e-3) / 20e-3
            i_hold = i_hold / 400e-12
            v_noise = v_noise / 5e-3
            i_noise = i_noise / 100e-12
            
            plt = self.get_channel_plot(chan)
            plt.setLabels(left=("Ch %d" % chan))
            for data,symbol in [(np.zeros_like(times), 'o'), (v_hold, 't'), (i_hold, 'x'), (v_noise, 't1'), (i_noise, 'x')]:
                brushes = np.where(np.abs(data) > 1.0, fail_brush, pass_brush)
                plt.plot(times, data, pen=None, symbol=symbol, symbolPen=None, symbolBrush=brushes)

        for i in recs.keys():
            start = (recs[i][0].start_time - start_time).seconds - 1
            stop = (recs[i][-1].start_time - start_time).seconds + 1
            pip_param = self.params.child('Pipette %d' % (i+1))
            pip_param.set_time_range(start, stop)
            
            got_data = len(recs[i]) > 2
            pip_param['got data'] = got_data
        
    def add_pipette(self, channel, status=None, **kwds):
        elec = PipetteParameter(self, channel, status=status, **kwds)
        self.params.addChild(elec, autoIncrementName=True)
        elec.child('channel').sigValueChanged.connect(self._pipette_channel_changed)
        elec.region.sigRegionChangeFinished.connect(self._pipette_region_changed)
        self._pipette_channel_changed(elec.child('channel'))
        
    def _pipette_channel_changed(self, param):
        plt = self.get_channel_plot(param.value())
        plt.addItem(param.parent().region)
        self._rename_pipettes()
        
    def _pipette_region_changed(self):
        self._rename_pipettes()
        
    def _rename_pipettes(self):
        # sort electrodes by channel
        elecs = {}
        for elec in self.params.children():
            elecs.setdefault(elec['channel'], []).append(elec)
        
        for chan in elecs:
            # sort all electrodes on this channel by start time
            chan_elecs = sorted(elecs[chan], key=lambda e: e.region.getRegion()[0])
            
            # assign names
            for i,elec in enumerate(chan_elecs):
                # rename all first to avoid name colisions
                elec.setName('rename%d' % i)
            for i,elec in enumerate(chan_elecs):
                # If there are multiple electrodes on this channel, then 
                # each extra electrode increments its name by the number of 
                # headstages (for example, on AD channel 3, the first electrode
                # is called "Electrode 4", and on an 8-headstage system, the
                # second electrode will be "Electrode 12").
                e_id = (chan+1) + (i*config.n_headstages)
                elec.id = e_id
                elec.setName('Pipette %d' % e_id)

    def save(self):
        state = {}
        for elec in self.params.children():
            rgn = elec.region.getRegion()
            if self.start_time is None:
                start = None
                stop = None
            else:
                start = self.start_time + datetime.timedelta(seconds=rgn[0])
                stop = self.start_time + datetime.timedelta(seconds=rgn[1])
            
            state[elec.id] = OrderedDict([
                ('pipette_status', elec['status']),
                ('got_data', elec['got data']),
                ('ad_channel', elec['channel']),
                ('patch_start', start),
                ('patch_stop', stop),
                ('cell_labels', {'biocytin': '', 'red': '', 'green': '', 'blue': ''}),
                #('cell_qc', {'holding': None, 'access': None, 'spiking': None}),
                ('target_layer', elec['target layer']),
                ('morphology', elec['morphology']),
                ('internal_solution', elec['internal']),
                ('internal_dye', elec['internal dye']),
                ('synapse_to', None),
                ('gap_to', None),
                ('notes', ''),
            ])
        return state


class PipetteParameter(pg.parametertree.parameterTypes.GroupParameter):
    def __init__(self, ui, channel, start=None, stop=None, status=None, got_data=None, internal=None, internal_dye=None, target_layer=None):
        self.ui = ui
        params = [
            {'name': 'channel', 'type': 'list', 'values': ui.list_channels()},
            {'name': 'status', 'type': 'list', 'values': ['No seal', 'Low seal', 'GOhm seal', 'Technical failure', 'No attempt', 'Not recorded']},
            {'name': 'got data', 'type': 'bool'},
            {'name': 'internal', 'type': 'list', 'values': [''] + constants.INTERNAL_RECIPES},
            {'name': 'internal dye', 'type': 'list', 'values': [''] + constants.INTERNAL_DYES},
            {'name': 'target layer', 'type': 'list', 'values': [''] + constants.LAYERS},
            {'name': 'morphology', 'type': 'list', 'values': {'':'', 'pyramidal':'pyr', 'non-pyramidal':'nonpyr', 'ambiguous':'?', 'no morphology':'x'}},
        ]
        pg.parametertree.parameterTypes.GroupParameter.__init__(self, name="Pipette?", children=params, removable=True)
        self.child('got data').sigValueChanged.connect(self._got_data_changed)
        
        region = [0, 1]
        if start is not None:
            region[0] = start
        if stop is not None:
            region[1] = stop

        self.region = pg.LinearRegionItem(region)
        self.region.setZValue(-10)
        
        if channel is not None:
            self['channel'] = channel
        if status is not None:
            self['status'] = status
        if got_data is not None:
            self['got data'] = got_data
        if internal is not None:
            self['internal'] = internal
        if internal_dye is not None:
            self['internal dye'] = internal_dye
        if target_layer is not None:
            self['target layer'] = target_layer
    
    def _got_data_changed(self):
        if self['got data'] is True:
            self.region.setBrush((0, 255, 0, 100))
        else:
            self.region.setBrush((0, 0, 255, 100))
        self.region.update()  # pg bug

    def set_time_range(self, start, stop):
        self.region.setRegion([start, stop])
