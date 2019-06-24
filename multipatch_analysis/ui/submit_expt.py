"""
Script used to submit completed experiment to database.
"""

from datetime import datetime

import acq4.util.Canvas, acq4.util.DataManager
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

from .. import metadata_submission, config
from . import ui


class ExperimentSubmitUi(QtGui.QWidget):
    def __init__(self):
        self.path = None
        
        QtGui.QWidget.__init__(self)
        self.setWindowTitle("Experiment metadata QC")
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.hsplit = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.hsplit, 0, 0)

        self.left_panel = QtGui.QWidget()
        self.left_layout = QtGui.QVBoxLayout()
        self.left_layout.setContentsMargins(0, 0, 0, 0)
        self.left_panel.setLayout(self.left_layout)
        self.hsplit.addWidget(self.left_panel)

        self.file_tree = FileTreeWidget(self)
        self.file_tree.itemSelectionChanged.connect(self.selection_changed)
        self.file_tree.itemDoubleClicked.connect(self.load_clicked)
        self.left_layout.addWidget(self.file_tree)
        
        self.ctrl_widget = QtGui.QWidget()
        self.ctrl_layout = QtGui.QGridLayout()
        self.ctrl_widget.setLayout(self.ctrl_layout)
        self.ctrl_layout.setContentsMargins(0, 0, 0, 0)
        self.left_layout.addWidget(self.ctrl_widget)
        
        row = self.ctrl_layout.rowCount()
        self.load_btn = QtGui.QPushButton('load files')
        self.load_btn.clicked.connect(self.load_clicked)
        self.ctrl_layout.addWidget(self.load_btn, row, 0)
        
        self.submit_btn = QtGui.QPushButton('submit...')
        self.submit_btn.clicked.connect(self.submit_clicked)
        self.submit_btn.setEnabled(False)
        self.ctrl_layout.addWidget(self.submit_btn, row, 1)

        self.vsplit = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.hsplit.addWidget(self.vsplit)
        self.hsplit.setSizes([600, 700])
        
        self.canvas = acq4.util.Canvas.Canvas(allowTransforms=False)
        self.vsplit.addWidget(self.canvas)
        
        self.timeline = ui.ExperimentTimeline()
        self.vsplit.addWidget(self.timeline)
        self.vsplit.setSizes([300, 600])
        
        self.submit_window = SubmitWindow()
        
    def set_path(self, path):
        self.path = path
        if path.info().get('dirType', None) != 'Site':
            raise Exception("Requested path is not a site directory: %s" % (path.name()))
        self.file_tree.set_path(path)
        self.timeline.load_site(path)
        self.canvas.clear()

    def load_clicked(self):
        sel = self.file_tree.selectedItems()
        for item in sel:
            fh = item.fh
            if isinstance(item, NwbTreeItem):
                self.timeline.load_nwb(fh)
            else:
                self.canvas.addFile(fh)

    def submit_clicked(self):
        sel = self.file_tree.selectedItems()[0]
        sub = sel.submission()
        self.submit_window.update_info(sub)
        self.submit_window.show()
        self.submit_window.activateWindow()        

    def selection_changed(self):
        sel = self.file_tree.selectedItems()
        sub = len(sel) == 1 and sel[0].is_submittable
        self.submit_btn.setEnabled(sub)


class FileTreeWidget(pg.TreeWidget):
    def __init__(self, ui):
        pg.TreeWidget.__init__(self)
        self.ui = ui
        self.path = None
        self.setColumnCount(3)
        self.setHeaderLabels(['file', 'category', 'metadata'])
        self.setSelectionMode(self.ExtendedSelection)
        self.setDragDropMode(self.NoDragDrop)
        
        # attempts to retain background colors on selected items:
        #self.setAllColumnsShowFocus(False)
        #self.itemSelectionChanged.connect(self._selection_changed)
        #self.style_delegate = StyleDelegate(self)
        #self.setItemDelegateForColumn(1, self.style_delegate)

    def set_path(self, path):
        self.path = path
        self._reload_file_tree()
        
    def _reload_file_tree(self):
        self.clear()
        
        dh = self.path.parent()
        root = self.invisibleRootItem()
        self._fill_tree(dh, root)
        
    def _fill_tree(self, dh, root):
        self.items = {}
        for fname in dh.ls():
            fh = dh[fname]
            if fh.isDir() and fh is not self.path:
                # exclude everything outside the selected site
                continue
            item = self._make_item(fh)
            self.items[fh] = item
            if hasattr(item, 'type_selected'):
                item.type_selected.connect(self._item_type_selected)
            
            root.addChild(item)
            item.setExpanded(True)
            item.fh = fh
            if fh.isDir():
                self._fill_tree(fh, item)
        
        for i in range(3):
            self.resizeColumnToContents(i)
        
    def _make_item(self, fh):
        info = fh.info()
        objtyp = info.get('__object_type__')
        
        if fh.isDir():
            dirtyp = info.get('dirType', None)
            dtyps = {'Experiment': ExperimentTreeItem, 'Slice': SliceTreeItem, 'Site': SiteTreeItem}
            if dirtyp in dtyps:
                return dtyps[dirtyp](self.ui, fh)
        
        if objtyp in ['ImageFile', 'MetaArray']:
            return ImageTreeItem(self.ui, fh)
        
        elif fh.shortName().lower().endswith('.nwb'):
            return NwbTreeItem(self.ui, fh)
        
        elif fh.shortName().startswith('MultiPatch_'):
            return MPLogTreeItem(self.ui, fh)
        
        item = TypeSelectItem(self.ui, fh, ['ignore'], 'ignore')
        return item

    def _item_type_selected(self, item, typ):
        for item in self.selectedItems():
            item.set_type(typ)
        self.resizeColumnToContents(1)

    ###### attempts to retain background colors on selected items:
    #def _selection_changed(self):
        ## Only select first column
        #try:
            #self.blockSignals(True)
            #for i in self.selectionModel().selectedIndexes():
                #if i.column() != 0:
                    #self.selectionModel().select(i, QtGui.QItemSelectionModel.Deselect)
        #finally:
            #self.blockSignals(False)

    #def mousePressEvent(self, ev):
        #if ev.button() == QtCore.Qt.RightButton:
            #print('press')
            #ev.accept()
        #else:
            #pg.TreeWidget.mousePressEvent(self, ev)

    #def mouseReleaseEvent(self, ev):
        #if ev.button() == QtCore.Qt.RightButton:
            #index = self.indexAt(ev.pos())
            #item, col = self.itemFromIndex(index)
            #print('release', item, col)
            #self._itemClicked(item, col)
        #else:
            #pg.TreeWidget.mouseReleaseEvent(self, ev)


#class StyleDelegate(QtGui.QStyledItemDelegate):
    #def __init__(self, table):
        #QtGui.QStyledItemDelegate.__init__(self)
        #self.table = table
    
    #def paint(self, painter, option, index):
        ##print(index.row(), index.column())
        #QtGui.QStyledItemDelegate.paint(self, painter, option, index)


class ExperimentTreeItem(pg.TreeWidgetItem):
    def __init__(self, ui, fh):
        self.fh = fh
        pg.TreeWidgetItem.__init__(self, [fh.shortName()])


class SliceTreeItem(pg.TreeWidgetItem):
    def __init__(self, ui, fh):
        self.fh = fh
        self.is_submittable = False
        
        #in_db = database.slice_from_timestamp(datetime.fromtimestamp(fh.info()['__timestamp__']))
        #if len(in_db) == 0:
            #status = "NOT SUBMITTED"
        #else:
            #status = "submitted"
        pg.TreeWidgetItem.__init__(self, [fh.shortName(), ''])


class SiteTreeItem(pg.TreeWidgetItem):
    def __init__(self, ui, fh):
        self.fh = fh
        self.ui = ui
        self.is_submittable = True

        pg.TreeWidgetItem.__init__(self, [fh.shortName(), ''])
        
    def submission(self):
        pips = self.ui.timeline.save()
        files = self.list_files()
        
        return metadata_submission.ExperimentMetadataSubmission(
            site_dh=self.fh, 
            files=files,
            pipettes=pips,
        )
    
    def list_files(self):
        """Generate a structure describing all files associated with this 
        experiment (site) and its parent slice.
        """
        files = []
        slice_dir = self.fh.parent()
        slice_item = self.parent()
        if slice_item is None:
            slice_item = self.treeWidget().invisibleRootItem()
            
        for parent in [slice_item, self]:
            childs = [parent.child(i) for i in range(parent.childCount())]
            for item in childs:
                if not item.fh.isFile():
                    continue
                typ = item.type()
                if typ == 'ignore':
                    continue
                files.append({'path': item.fh.name(relativeTo=slice_dir), 'category': typ})
        return files
        


class TypeSelectItem(pg.TreeWidgetItem):
    """TreeWidgetItem with a type selection menu in the second column.
    """
    class Signals(QtCore.QObject):
        type_selected = QtCore.Signal(object, object)
    
    def __init__(self, ui, fh, types, current_type):
        self.is_submittable = False
        self.fh = fh
        self._sigprox = ImageTreeItem.Signals()
        self.type_selected = self._sigprox.type_selected
        self.types = types
        pg.TreeWidgetItem.__init__(self, [fh.shortName(), '', ''])

        self.menu = QtGui.QMenu()
        for typ in self.types:
            act = self.menu.addAction(typ, self._type_selected)
        
        self.set_type(current_type)

    def _type_selected(self):
        action = self.treeWidget().sender()
        text = str(action.text()).strip()
        self.set_type(text)
        self.type_selected.emit(self, text)
            
    def set_type(self, typ):
        self.setText(1, typ)
        if typ == 'ignore':
            self.setBackground(1, pg.mkColor(0.9))
        else:
            self.setBackground(1, pg.mkColor('w'))

    def type(self):
        return self.text(1)

    def itemClicked(self, col):
        if col != 1:
            return
        tw = self.treeWidget()
        x = tw.header().sectionPosition(col)
        y = tw.header().height() + tw.visualItemRect(self).bottom()
        self.menu.popup(tw.mapToGlobal(QtCore.QPoint(x, y)))
        return None
        

class NwbTreeItem(TypeSelectItem):
    def __init__(self, ui, fh):
        types = ['ignore', 'MIES physiology']
        if fh.parent().info().get('dirType') == 'Site':
            typ = 'MIES physiology'
        else:
            typ = 'ignore'
        TypeSelectItem.__init__(self, ui, fh, types, typ)
    

class MPLogTreeItem(TypeSelectItem):
    def __init__(self, ui, fh):
        types = ['ignore', 'Multipatch log']
        if fh.parent().info().get('dirType') == 'Site':
            typ = types[1]
        else:
            typ = 'ignore'
        TypeSelectItem.__init__(self, ui, fh, types, typ)
    

class ImageTreeItem(TypeSelectItem):
    def __init__(self, ui, fh):
        info = fh.info()
        obj = fh.info().get('objective', '')

        # Make initial guess on image type
        typ = 'ignore'
        ptype = fh.parent().info().get('dirType')
        if ptype == 'Site':
            typ = 'recording site'
        elif ptype == 'Slice':
            if obj.startswith('40x'):
                typ = 'slice quality stack'
            if obj.startswith('4x'):
                typ = 'slice anatomy'

        types = ['ignore', 'slice anatomy', 'slice quality stack', 'recording site']
        TypeSelectItem.__init__(self, ui, fh, types, typ)
        
        self.setText(2, obj)
        illumination = info.get('illumination', {})
        if illumination is not None:
            colors = list(illumination.keys())
        else:
            colors = []
        if len(colors) == 0:
            color = 'w'
        elif len(colors) > 1:
            color = 'y'
        else:
            color = {'infrared': (255, 230, 230), 'green': (230, 255, 230), 'blue': (230, 230, 255), 'uv': (255, 220, 255)}[colors[0]]
        self.setBackground(2, pg.mkColor(color))
            

class SubmitWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.resize(800, 800)
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.message_text = QtGui.QTextBrowser()
        self.message_text.setOpenExternalLinks(True)
        self.layout.addWidget(self.message_text, 0, 0)
        
        self.info_tree = pg.DataTreeWidget()
        self.info_tree.setColumnHidden(1, True)
        self.layout.addWidget(self.info_tree, 0, 1)
        
        self.submit_btn = QtGui.QPushButton('submit!')
        self.submit_btn.clicked.connect(self.submit)
        self.layout.addWidget(self.submit_btn, 1, 1)
        
    def update_info(self, submission):
        self.submission = submission
        errors, warnings = submission.check()
        messages = ''
        for name, msgs in [('errors', errors), ('warnings', warnings)]:
            if len(msgs) == 0:
                messages += "<h3>No %s.</h3><br><br>" % name
            else:
                messages += "<h3>%s:</h3></br>\n<ul>\n" % name.capitalize()
                for msg in msgs:
                    messages += "<li>" + msg + "<br>\n"
                messages += "</ul><br><br>\n"
        
        self.message_text.setHtml(messages)
        
        summary = submission.summary()
        self.info_tree.setData(summary)
        
    def submit(self):
        if len(self.submission.check()[0]) > 0:
            raise Exception("Can't submit; experiment has errors.")
        self.submission.submit()
        self.hide()
        
        


def submit(data):
    print "Submitting ", data
    session = Session()
    


        
if __name__ == '__main__':
    import sys
    app = pg.mkQApp()
    pg.dbg()
    
    path = acq4.util.DataManager.getDirHandle(sys.argv[1])
    ui = ExperimentSubmitUi()
    ui.resize(1600, 1000)
    ui.show()
    ui.set_path(path)
    
    if sys.flags.interactive == 0:
        app.exec_()
