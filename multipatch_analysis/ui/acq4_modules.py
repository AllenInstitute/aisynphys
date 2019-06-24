import acq4.pyqtgraph as pg
from acq4.pyqtgraph.Qt import QtCore, QtGui
from acq4.modules.Module import Module
from acq4.Manager import getManager
import acq4.util.DataManager as DataManager

# from acq4.analysis.modules.MosaicEditor import MosaicEditor
# from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.modules.MosaicEditor import MosaicEditor
from acq4.util.Canvas.items.MarkersCanvasItem import MarkersCanvasItem 
from acq4.util.Canvas.items.MultiPatchLogCanvasItem import MultiPatchLogCanvasItem 
from . import submit_expt
from . import multipatch_nwb_viewer
from .. import lims
from . import vimaging
import os
import shutil
import json
import urllib
from datetime import datetime
import time
import yaml

from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.util.Canvas.items import registerItemType
from . import dashboard


class MultipatchSubmissionModule(Module):
    """Allows multipatch data submission UI to be invoked as an ACQ4 module.
    
    This is primarily to ensure that metadata changes made via ACQ4 are immediately
    available to (and do noty collide with) the submission tool, and vice-versa. 
    """
    moduleDisplayName = "MP Submission Tool"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = submit_expt.ExperimentSubmitUi()
        self.ui.resize(1600, 900)
        self.ui.show()
        
        self.load_from_dm_btn = QtGui.QPushButton("load from data manager")
        self.ui.left_layout.insertWidget(0, self.load_from_dm_btn)
        self.load_from_dm_btn.clicked.connect(self.load_from_dm_clicked)
        
    def load_from_dm_clicked(self):
        man = getManager()
        sel_dir = man.currentFile
        self.ui.set_path(sel_dir)
        
    def window(self):
        return self.ui

        
class NWBViewerModule(Module):
    """ACQ module for browsing data in NWB files.
    """
    moduleDisplayName = "MP NWB Viewer"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = multipatch_nwb_viewer.MultipatchNwbViewer()
        self.ui.resize(1600, 900)
        self.ui.show()
        
        self.load_from_dm_btn = QtGui.QPushButton("load from data manager")
        self.ui.vsplit.insertWidget(0, self.load_from_dm_btn)
        self.load_from_dm_btn.clicked.connect(self.load_from_dm_clicked)
        
    def load_from_dm_clicked(self):
        man = getManager()
        filename = man.currentFile.name()
        nwb = self.ui.load_nwb(filename)
        
    def window(self):
        return self.ui


class MultiPatchMosaicEditorExtension(QtGui.QWidget):
    def __init__(self, mosaic_editor):
        self.mosaic_editor = mosaic_editor

        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.load_btn = QtGui.QPushButton("Download 20x")
        self.layout.addWidget(self.load_btn, 0, 0)
        self.load_btn.clicked.connect(self.load_clicked)

        self.submit_btn = pg.FeedbackButton("Submit cell positions to LIMS")
        self.layout.addWidget(self.submit_btn, 1, 0)
        self.submit_btn.clicked.connect(self.submit)

        self.map_no_go = QtGui.QCheckBox("Cell map no go")
        self.layout.addWidget(self.map_no_go)

    def load_clicked(self):
        """
        Checks that directory is pointed at a slice
        Check LIMS for 20x biocytin image and if 20x image is saved locally.
        If image not saved locally then copies file from LIMS.
        Enables button to create a JSON of the cell locations
        """
        if self.base_dir.info()['dirType'] != 'Slice':
            raise Exception('No Slice Selected')
        self.base_path = self.base_dir.name()

        aff_image_path = lims.specimen_20x_image(self.specimen_name, treatment='Biocytin')
        if aff_image_path is None:
            raise Exception("No 20x biocytin image for specimen %s" % self.specimen_name)

        if not os.path.exists(aff_image_path):
            raise Exception("No image found at path '%s'" % aff_image_path)
        
        aff_image_name = os.path.split(aff_image_path)[-1]
        save_path = os.path.join(self.base_path, aff_image_name)
            
        if not os.path.exists(save_path):
            with pg.BusyCursor():
                shutil.copy2(aff_image_path, save_path)
                self.base_dir.indexFile(aff_image_name)
        
    
    @property
    def base_dir(self):
        return self.mosaic_editor.ui.fileLoader.baseDir()
    
    @property
    def specimen_name(self):
        if self.base_dir.info()['dirType'] != 'Slice':
            raise Exception('No Slice Selected')
        return self.base_dir.info()['specimen_ID'].strip()

    def submit(self):
        try:
            self.create_lims_json()
            self.save_json_and_trigger()
            self.submit_btn.success()
        except:
            self.submit_btn.failure()
            raise

    def create_lims_json(self):
        """
        Checks that markers are loaded into the canvas
        Gets timestamp
        Creates dictionary of cell locations and enables the save button
        If the cluster is unmappable (no stained cells, damage, etc.) the json will submit cell_map_no_go to specimen tags
        """

        cell_map_no_go = self.map_no_go.isChecked()
        if cell_map_no_go is True:
            self.data = {}
            self.data['cell map no go'] = cell_map_no_go
        else:
            items = self.mosaic_editor.canvas.items
            image_20 = lims.specimen_20x_image(self.specimen_name, treatment='Biocytin')

            # Find item that marks cell positions
            markers = [i for i in items if isinstance(i, MarkersCanvasItem)]
            if len(markers) != 1:
                raise Exception("Must have exactly 1 Markers item in the canvas.")
            
            # Find the 20x image that was loaded in to the editor and make sure it's the right
            # one for this specimen
            image = None
            if image_20 is not None:
                for i in items:
                    try:
                        if os.path.split(i.opts['handle'].name())[-1] == os.path.split(image_20)[-1]:
                            image = i
                    except AttributeError:
                        pass
            # if image is None:
            #     raise Exception("Could not find 20x image loaded into mosaic editor.")

            # Infer the site folder to submit from based on the files that have been loaded.
            # This is really indirect, but it seems to work.
            sites = set()
            for item in items:
                handle = item.opts.get('handle')
                if handle is None:
                    continue
                parent = handle.parent()
                if parent.shortName().startswith('site_'):
                    sites.add(parent)
            if len(sites) == 0:
                raise Exception("No files loaded from site; can't determine which site to submit.")
            if len(sites) > 1:
                raise Exception("Files loaded from multiple sites; can't determine which site to submit.")
            self.site_dh = list(sites)[0]

            # Grab the timestamp for this site; we need it to look up the cell cluster from LIMS
            ts = self.site_dh.info()['__timestamp__']

            # Find the LIMS CellCluster ID for this experiment 
            cluster_ids = lims.expt_cluster_ids(self.specimen_name, ts)
            if len(cluster_ids) == 0:
                raise Exception("No CellCluster found for %s %s" % (self.specimen_name, ts))
            if len(cluster_ids) > 1:
                raise Exception("Multiple CellClusters found for %s %s" % (self.specimen_name, ts))
            self.cluster_id = cluster_ids[0]
            self.cluster_name = lims.specimen_name(self.cluster_id)

            # Load up the piettes.yml file so we know which cells are considered "ephys pass"
            pipette_path = os.path.join(self.site_dh.name(), 'pipettes.yml')
            if not os.path.exists(pipette_path):
                raise Exception("Could not find pipettes.yml")
            pipette_log = yaml.load(open(pipette_path))

            # grab info about cell positions from Markers
            pipettes = markers[0].saveState()['markers']

            # Construct JSON to send to LIMS.
            self.data = {}  
            self.data['cells'] = []
            for cell in pipettes:
                cell_name = int(cell[0].split('_')[-1])
                p = cell[1]
                if image is None:
                    x_float = 0
                    y_float = 0
                else:
                    x_float = markers[0].graphicsItem().mapToItem(image.graphicsItem(), pg.Point(*p[:2])).x()
                    y_float = markers[0].graphicsItem().mapToItem(image.graphicsItem(), pg.Point(*p[:2])).y()
                self.data['cells'].append({      
                    'ephys_cell_id': cell_name,                  
                    'coordinates_20x': 
                    {
                    "x": int(round(x_float)),   # round the float and change to integer for technology requirements
                    "y": int(round(y_float))    # round the float and change to integer for technology requirements
                    },      
                    'ephys_qc_result': ('pass' if pipette_log[cell_name]['got_data'] == True else 'fail'),                 
                    'start_time_sec': time.mktime(pipette_log[cell_name]['patch_start'].timetuple())
                })

            # check that position values written to self.data are non-zero, if they are ask the user if it's ok to submit
            x_pos = all([cell['coordinates_20x']['x'] for cell in self.data['cells']]) != 0
            y_pos = all([cell['coordinates_20x']['y'] for cell in self.data['cells']]) != 0
            if x_pos is False or y_pos is False:
                ret = QtGui.QMessageBox.question(self, "Position Check", "x and/or y position values are 0.\n"
                    "Do you still want to submit?", QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel)
                # postionValueMsg = QMessageBox()
                # postionValueMsg.setText('x and/or y position values are 0')
                # positionValueMsg.setInformativeText('Do you still want to submit?')
                # positionValueMsg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                # answer = positionValueMsg.exec_()

                if ret == QtGui.QMessageBox.Cancel:
                    raise Exception ('Submission Cancelled')
        
        
    def save_json_and_trigger(self):
        """
        creates the json file and trigger file and moves them to the proper
        incoming folders for LIMS
        disables create json button
        """

        json_save_file = 'synphys_' + self.cluster_name + "_ephys_cell_cluster.json"

        # Ask LIMS where to put incoming files
        incoming_path = lims.get_incoming_dir(self.specimen_name)
        if incoming_path is None:
            raise Exception("Could not find LIMS incoming directory for specimen '%s'" % self.specimen_name)

        incoming_save_path = '//' + incoming_path.lstrip('/')
        if not os.path.isdir(incoming_save_path):
            raise Exception("LIMS incoming directory for specimen '%s' does not exist: '%s'" % (self.specimen_name, incoming_save_path))
        
        # Save payload json file
        incoming_json_save_path = incoming_save_path + '/' + json_save_file # using '/' on windows because it is a network path
        with open(incoming_json_save_path, 'w') as outfile:  
            json.dump(self.data, outfile)

        trigger_file = 'synphys_' + self.cluster_name + "_ephys_cell_cluster.mpc"

        trigger_path = lims.get_trigger_dir(self.specimen_name)
        if trigger_path is None:
            raise Exception("Could not find LIMS trigger directory for specimen '%s'" % self.specimen_name)
        trigger_path = '//' + trigger_path.lstrip('/')
        
        trigger_save_path = os.path.join(trigger_path, trigger_file)
        print trigger_save_path
        print incoming_json_save_path
        with open(trigger_save_path, 'w') as the_file:
            the_file.write("specimen_id: {}\n".format(self.cluster_id))     
            the_file.write("cells_info: '{}'\n".format(incoming_json_save_path))
        

MosaicEditor.addExtension("Multi Patch", {
    'type': 'ctrl',
    'builder': MultiPatchMosaicEditorExtension,
    'pos': ('top', 'Canvas'),
    'size': (600, 200),
})


class DashboardModule(Module):
    """ACQ module for monitoring pipeline status
    """
    moduleDisplayName = "MP Dashboard"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = dashboard.Dashboard(limit=config.get('limit', None), filter_defaults=config.get('filters', None))
        self.ui.resize(1600, 900)
        self.ui.show()
        
    def window(self):
        return self.ui


class VoltageImagingAnalysisModule(Module):
    """
    """
    moduleDisplayName = "Voltage Imaging"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = vimaging.VImagingAnalyzer()
        self.ui.resize(1600, 900)
        self.ui.show()
        
        self.load_from_dm_btn = QtGui.QPushButton("load from data manager")
        self.ui.layout.addWidget(self.load_from_dm_btn, self.ui.layout.rowCount(), 0)
        self.load_from_dm_btn.clicked.connect(self.load_from_dm_clicked)
        
    def load_from_dm_clicked(self):
        man = getManager()
        sel_dir = man.currentFile
        self.ui.load_data(sel_dir)
        
    def window(self):
        return self.ui


class VoltageImagingAnalysis2Module(Module):
    """
    """
    moduleDisplayName = "Voltage Imaging 2"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = vimaging.VImagingAnalyzer2()
        self.ui.resize(1600, 900)
        self.ui.show()
        
        self.load_from_dm_btn = QtGui.QPushButton("load from data manager")
        self.load_from_dm_btn.setParent(self.ui)
        self.load_from_dm_btn.resize(160, 30)
        self.load_from_dm_btn.show()
        self.load_from_dm_btn.clicked.connect(self.load_from_dm_clicked)
        
    def load_from_dm_clicked(self):
        man = getManager()
        sel_dir = man.currentFile
        self.ui.load_data(sel_dir)
        
    def window(self):
        return self.ui


class AffPyramidCanvasItem(CanvasItem):
    """For displaying AFF image pyramids
    """
    _typeName = "AFF Image Pyramid"
    
    def __init__(self, handle, **kwds):
        from affpyramid.ui import AffImageItem
        kwds.pop('viewRect', None)
        self.affitem = AffImageItem(handle.name())
        opts = {'movable': True, 'rotatable': True, 'handle': handle}
        opts.update(kwds)
        if opts.get('name') is None:
            opts['name'] = handle.shortName()
        opts['defaultUserTransform'] = {'scale': (0.36e-6, 0.36e-6)}            
        CanvasItem.__init__(self, self.affitem, **opts)

    @classmethod
    def checkFile(cls, fh):
        name = fh.shortName()
        if name.endswith('.aff'):
            return 10
        else:
            return 0

registerItemType(AffPyramidCanvasItem)
