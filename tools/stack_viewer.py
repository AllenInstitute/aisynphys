import os, sys, argparse
import pyqtgraph as pg
import numpy as np
from multipatch_analysis import lims
from affpyramid.ui import AffImageItem

# get specimen ID from command line
parser = argparse.ArgumentParser("Viewer for image stacks served from LIMS")
parser.add_argument("specimen", type=str, help="Specimen ID or name")
args = parser.parse_args()
try:
    args.specimen = int(args.specimen)
except ValueError:
    pass

# look up image file names
jp2_files = sorted(lims.specimen_images(args.specimen)[0]['file'], key=lambda f: os.path.split(f)[-1])
if sys.platform == 'win32':
    jp2_files = ['\\' + '\\'.join(f.split('/')) for f in jp2_files]
aff_files = [os.path.splitext(f)[0] + '.aff' for f in jp2_files]


# make a window with zoomable view box
app = pg.mkQApp()
window = pg.GraphicsLayoutWidget()
window.show()
view = window.addViewBox()

# view should have square pixels and y-axis pointig downward
view.setAspectLocked(True)
view.invertY()


# for displaying a particular z slice
aff_image_item = None
z_index = None
def set_z_index(zi):
    global aff_image_item, z_index
    zi = np.clip(zi, 0, len(aff_files)-1)
    z_index = zi
    if aff_image_item is not None:
        view.removeItem(aff_image_item)
    
    aff_image_item = AffImageItem(aff_files[int(zi)])
    view.addItem(aff_image_item)


set_z_index(len(aff_files) // 2)


# catch pgup/pgdn to set z slice

class KeyHandler(pg.QtCore.QObject):
    def __init__(self):
        pg.QtCore.QObject.__init__(self)
        self.pressed_keys = set()
        self.modifiers = 0
        self.timer = pg.QtCore.QTimer()
        self.timer.timeout.connect(self.update)
        self.timer.start(30)

    def eventFilter(self, obj, event):
        if event.type() not in (event.KeyPress, event.KeyRelease):
            return False

        if event.key() not in (pg.QtCore.Qt.Key_PageUp, pg.QtCore.Qt.Key_PageDown):
            return False

        if event.isAutoRepeat():
            return True

        if event.type() == event.KeyPress:
            self.pressed_keys.add(event.key())
        if event.type() == event.KeyRelease:
            self.pressed_keys.remove(event.key())
        self.modifiers = event.modifiers()

        return True

    def update(self):
        # Handle key presses on a timer so we get smooth scrolling across z planes

        zi = z_index

        dz = 1.0
        if int(self.modifiers & pg.QtCore.Qt.ControlModifier) > 0:
            dz /= 5.0
        if int(self.modifiers & pg.QtCore.Qt.ShiftModifier) > 0:
            dz *= 5.0

        if pg.QtCore.Qt.Key_PageUp in self.pressed_keys:
            zi -= dz
        if pg.QtCore.Qt.Key_PageDown in self.pressed_keys:
            zi += dz

        if zi == z_index:
            return

        set_z_index(zi)
        app.processEvents()


key_handler = KeyHandler()
window.installEventFilter(key_handler)


# Start event loop if the script is not running interactively
if sys.flags.interactive == 0:
    app.exec_()

