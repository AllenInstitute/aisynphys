from acq4.Manager import getManager
from acq4.pyqtgraph.Qt import QtGui, QtCore
import acq4.pyqtgraph as pg
import numpy as np
import scipy.ndimage as ndimage
import scipy.stats as stats


class VImagingAnalyzer(QtGui.QSplitter):
    def __init__(self):
        QtGui.QSplitter.__init__(self, QtCore.Qt.Horizontal)
        self.resize(800, 1000)
            
        self.params = pg.parametertree.Parameter(name='params', type='group', children=[
            dict(name='sequence', type='group', children=[
                dict(name='analysis', type='list', values=['dF / F', 'SNR', 'noise']),
            ]),
        ])
        self.ptree = pg.parametertree.ParameterTree()
        self.ptree.setParameters(self.params)
        self.params.child('sequence').sigTreeStateChanged.connect(self.update_sequence_analysis)        
        
        self.leftPanel = QtGui.QWidget()
        self.addWidget(self.leftPanel)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.leftPanel.setLayout(self.layout)
        
        self.layout.addWidget(self.ptree, 0, 0)

        # adds a gradient widget to the left panel; later linked to vb2
        self.grw = pg.GradientWidget()
        self.grw.loadPreset('bipolar')
        self.layout.addWidget(self.grw)

        self.gw = pg.GraphicsLayoutWidget()
        self.addWidget(self.gw)
        
        self.vb1 = self.gw.addViewBox()
        self.img1 = pg.ImageItem()
        self.vb1.addItem(self.img1)
        
        self.vb2 = self.gw.addViewBox(row=1, col=0)
        self.img2 = pg.ImageItem()
        self.vb2.addItem(self.img2)
        # links axes of viewboxes (fluor and dF)
        self.vb2.setXLink(self.vb1)
        self.vb2.setYLink(self.vb1)
        
        for vb in (self.vb1, self.vb2):
            vb.invertY()
            vb.setAspectLocked(True)

        self.plt1 = self.gw.addPlot(row=2, col=0)
        self.plt1.setLabels(bottom=('time', 's'))
        self.plt1_items = []
        
        self.rois = [pg.EllipseROI([0, 0], [10, 10], pen=(i,4)) for i in range(2)]
        for roi in self.rois:
            self.vb1.addItem(roi)
            roi.sigRegionChangeFinished.connect(self.roi_changed)
        self.ignore_roi_change = False
            
        self.plt2 = self.gw.addPlot(row=3, col=0)
        self.plt2.setXLink(self.plt1)
        self.plt2.setLabels(bottom=('time', 's'))
            
        self.noise_time_rgn = pg.LinearRegionItem([0.01, 0.02], brush=(255, 0, 0, 30))
        self.base_time_rgn = pg.LinearRegionItem([0.025, 0.95], brush=(0, 255, 0, 30))
        self.test_time_rgn = pg.LinearRegionItem([0.1, 0.11], brush=(0, 0, 255, 30))
        for time_rgn in (self.base_time_rgn, self.test_time_rgn, self.noise_time_rgn):
            self.plt1.addItem(time_rgn)
            time_rgn.sigRegionChangeFinished.connect(self.time_rgn_changed)
            time_rgn.sigRegionChangeFinished.connect(self.update_sequence_analysis)
        self.clamp_plots = []
        
        self.plt3 = self.gw.addPlot(row=4, col=0)
        self.plt3.setLabels(left="dF / F", bottom="trial")
        
        self.show()

    def load_data(self, seqDir):
        with pg.BusyCursor():
            self._load_data(seqDir)

    def _load_data(self, seqDir):
        man = getManager()
        model = man.dataModel

        # read all image data
        self.img_data = model.buildSequenceArray(
            seqDir, 
            lambda dh: dh['Camera']['frames.ma'].read(),
            join=False).asarray()
        seqParams = list(model.listSequenceParams(seqDir).items())
        if self.img_data.ndim == 1:
            self.img_data = self.img_data[np.newaxis, :]
            seqParams.insert(0, (None, [0]))

        transpose = seqParams[0][0] == ('protocol', 'repetitions')
        if transpose:
            self.img_data = np.swapaxes(self.img_data, 0, 1)
            seqParams = seqParams[::-1]
        self.seqParams = seqParams

        if seqParams[0][0] is None:
            self.seqColors = [pg.mkColor('w')]
        else:
            nSeq = len(seqParams[0][1])
            if nSeq == 2:
                # Special case: add in difference between two sequence trials
                seqParams[0] = (seqParams[0][0], list(seqParams[0][1]) + [np.mean(seqParams[0][1])])
                nSeq = 3
                img_data = np.empty((3,) + self.img_data.shape[1:], dtype=self.img_data.dtype)
                img_data[:2] = self.img_data
                for i in range(img_data.shape[1]):
                    minlen = min(self.img_data[0,i].shape[0], self.img_data[1,i].shape[0])
                    img_data[2,i] = self.img_data[0,i][:minlen].copy()
                    img_data[2,i]._data = self.img_data[0,i][:minlen].asarray().astype('float') - self.img_data[1,i][:minlen].asarray()
                self.img_data = img_data
                
                
            self.seqColors = [pg.intColor(i, nSeq*1.6) for i in range(nSeq)]

        # cull out truncated recordings :(
        self.img_data = [[d['Time':0:500e-3] for d in row if d.xvals('Time')[-1] > 150e-3] for row in self.img_data]

        # crop / concatenate
        img_len = min([min([d.shape[0] for d in row]) for row in self.img_data])
        self.img_data = [[d[:img_len] for d in row] for row in self.img_data]
        self.img_arrays = [np.concatenate([d.asarray()[np.newaxis, ...] for d in row], axis=0) for row in self.img_data]

        # average (across trials???)
        self.img_mean = [img_arr.mean(axis=0) for img_arr in self.img_arrays]

        for p in self.clamp_plots:
            self.plt2.removeItem(p)
            
        # read all clamp data
        first_subdir = seqDir[seqDir.ls()[0]]
        clamp_name = 'Clamp1'
        if not first_subdir[clamp_name + '.ma'].exists():
            clamp_name = 'Clamp2'
        
        clamp_file = first_subdir[clamp_name + '.ma']
        if clamp_file.exists():
            self.clamp_mode = model.getClampMode(clamp_file)
            chan = 'command' if self.clamp_mode == 'VC' else 'primary'
            self.clamp_data = model.buildSequenceArray(
                seqDir, 
                lambda dh: dh[clamp_name + '.ma'].read()['Channel': chan],
                join=False).asarray()
            if self.clamp_data.ndim == 1:
                self.clamp_data = self.clamp_data[np.newaxis, :]
            if transpose:
                self.clamp_data = np.swapaxes(self.clamp_data, 0, 1)

            self.plt2.setLabels(left=('Vm', 'V'))
            
            for i in range(self.clamp_data.shape[0]):
                for j in range(self.clamp_data.shape[1]):
                    trace = self.clamp_data[i,j]
                    pen = self.seqColors[i]
                    p = self.plt2.plot(trace.xvals('Time'), trace.asarray(), antialias=True, pen=pen)
                    self.clamp_plots.append(p)
            self.plt2.show()
        else:
            self.clamp_mode = None
            self.plt2.hide()

        self.img_t = self.img_data[0][0].xvals('Time')

        self.img1.setImage(self.img_mean[-1].mean(axis=0))

        self.roi_changed(None)
        self.time_rgn_changed(None)
        self.update_sequence_analysis()

    def roi_changed(self, roi):
        if self.ignore_roi_change:
            return
        
        # Update other ROI to match
        roi1, roi2 = self.rois
        try:
            if roi is not None:
                self.ignore_roi_change = True
                other_roi = roi2 if roi is roi1 else roi1
                pos1 = roi.pos()
                pos2 = other_roi.pos()
                pos2.setY(pos1.y())
                other_roi.setPos(pos2)
                other_roi.setSize(roi.size())            
        finally:
            self.ignore_roi_change = False
        
        
        for item in self.plt1_items:
            self.plt1.removeItem(item)

        
        for i, img_data in enumerate(self.img_arrays):
            color = self.seqColors[i]
            color2 = pg.mkColor(color)
            color2.setAlpha(40)

            rgn1 = roi1.getArrayRegion(img_data, self.img1, axes=(2, 3)).mean(axis=2).mean(axis=2)
            rgn2 = roi2.getArrayRegion(img_data, self.img1, axes=(2, 3)).mean(axis=2).mean(axis=2)
            dif = rgn1 - rgn2
            
            difmean = dif.mean(axis=0)
            baseline = np.median(difmean[:10])
            # plot individual examples only for last parameter
            if i == len(self.img_arrays)-1:
                for j in range(dif.shape[0]):
                    offset = baseline - np.median(dif[j, :10])
                    self.plt1_items.append(self.plt1.plot(self.img_t, dif[j] + offset, pen=color2, antialias=True))
            # plot average
            self.plt1_items.append(self.plt1.plot(self.img_t, difmean, pen=color, antialias=True))
            self.plt1_items[-1].setZValue(10)
            
        self.update_sequence_analysis()
        
    def time_rgn_changed(self, rgn):
        # make sure noise and test regions have the same width
        nr = self.noise_time_rgn.getRegion()
        tr = self.test_time_rgn.getRegion()
        if rgn is self.test_time_rgn:
            self.noise_time_rgn.setRegion([nr[0], nr[0] + tr[1]-tr[0]])
        elif rgn is self.noise_time_rgn:
            self.test_time_rgn.setRegion([tr[0], tr[0] + nr[1]-nr[0]])
            
        img_data = self.img_arrays[-1]
        img_t = self.img_data[-1][0].xvals('Time')
        
        base_starti, base_stopi, test_starti, test_stopi = self.time_indices(img_t)[:4]

        base = img_data[:, base_starti:base_stopi].mean(axis=0).mean(axis=0)
        test = img_data[:, test_starti:test_stopi].mean(axis=0).mean(axis=0)

        dff = ndimage.median_filter(test - base, 3) # original median radius was 10
        self.img2.setImage(dff) # mean across window, then mean across trials, then subtract
        self.img2.setLookupTable(self.grw.getLookupTable(65535)) # set LUT to be the gradient in the left panel widget

    def time_indices(self, time_vals):
        base_start, base_stop = self.base_time_rgn.getRegion()
        base_starti = np.argwhere(time_vals >= base_start)[0,0]
        base_stopi = np.argwhere(time_vals >= base_stop)[0,0]
        
        test_start, test_stop = self.test_time_rgn.getRegion()
        test_starti = np.argwhere(time_vals >= test_start)[0,0]
        test_stopi = np.argwhere(time_vals >= test_stop)[0,0]
        
        noise_start, noise_stop = self.noise_time_rgn.getRegion()
        noise_starti = np.argwhere(time_vals >= noise_start)[0,0]
        noise_stopi = np.argwhere(time_vals >= noise_stop)[0,0]
        
        return base_starti, base_stopi, test_starti, test_stopi, noise_starti, noise_stopi

    def update_sequence_analysis(self):
        analysis = self.params['sequence', 'analysis']
        xvals = []
        yvals = []
        noise = []
        brushes = []
        avg_x = []
        avg_y = []
        avg_noise = []

        for i in range(len(self.img_data)):

            img_data = self.img_arrays[i]
            img_t = self.img_data[i][0].xvals('Time')
            base_starti, base_stopi, test_starti, test_stopi, noise_starti, noise_stopi = self.time_indices(img_t)
            dff = self.measure_dff(img_data, img_t, (base_starti, base_stopi, test_starti, test_stopi))
            
            x = self.seqParams[0][1][i]
            xvals.extend([x] * len(dff))
            
            yvals.extend(list(dff))
            
            if analysis in ('SNR', 'noise'):
                noise_dff = self.measure_dff(img_data, img_t, (base_starti, base_stopi, noise_starti, noise_stopi))
                noise.extend(list(noise_dff))
                avg_noise.append(noise_dff.mean())
            
            color = pg.mkColor(self.seqColors[i])
            color.setAlpha(50)
            brushes.extend([pg.mkBrush(color)] * len(dff))
            avg_x.append(x)
            avg_y.append(dff.mean())

        if analysis == 'SNR':
            noise = np.std(noise)
            yvals = np.array(yvals) / noise
            avg_y = np.array(avg_y) / noise
        elif analysis == 'noise':
            yvals = noise
            avg_y = avg_noise


        self.plt3.clear()
        self.plt3.setLabels(left=analysis)
        if self.seqParams[0][0] is None:
            self.plt3.plot(range(len(xvals)), yvals, pen='w')
        else:
            self.plt3.plot(xvals, yvals, symbol='o', pen=None, symbolBrush=brushes, symbolPen=None)
            self.plt3.plot(avg_x, avg_y, pen='w', antialias=True)
            
            lin = stats.linregress(avg_x, avg_y)
            self.plt3.setTitle("slope: %0.2g" % lin[0])
    
    def measure_dff(self, img_data, img_t, time_indices):
        base_starti, base_stopi, test_starti, test_stopi = time_indices
        
        base = img_data[:, base_starti:base_stopi]
        base_mean = base.mean(axis=1)
        test = img_data[:, test_starti:test_stopi]
        test_mean = test.mean(axis=1)
        
        roi1, roi2 = self.rois   # roi1 is signal, roi2 is background
        base_rgn1 = roi1.getArrayRegion(base_mean, self.img1, axes=(1, 2))
        base_rgn2 = roi2.getArrayRegion(base_mean, self.img1, axes=(1, 2))
        test_rgn1 = roi1.getArrayRegion(test_mean, self.img1, axes=(1, 2))
        test_rgn2 = roi2.getArrayRegion(test_mean, self.img1, axes=(1, 2))

        # Use the temporal profile in roi2 in order to remove changes in LED brightness over time
        # Then use the difference between baseline and test time regions to determine change in fluorescence
        baseline1 = base_rgn1.mean(axis=1).mean(axis=1)
        signal1 = test_rgn1.mean(axis=1).mean(axis=1)
        baseline2 = base_rgn2.mean(axis=1).mean(axis=1)
        signal2 = test_rgn2.mean(axis=1).mean(axis=1)
        dff = ((signal1-signal2) - (baseline1-baseline2)) / (baseline1-baseline2)

        return dff
    
    def closeEvent(self, ev):
        self.img_data = None
        self.img_arrays = None
        self.img_mean = None
        self.clamp_data = None
        self.clamp_mode = None


class VImagingAnalyzer2(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)

        self.layout = QtGui.QVBoxLayout()
        self.setLayout(self.layout)
        self.imv1 = pg.ImageView()
        self.layout.addWidget(self.imv1)

        self.plt1 = pg.PlotWidget()
        self.layout.addWidget(self.plt1)

        self.base_time_rgn = pg.LinearRegionItem([0.07, 0.099])
        self.test_time_rgn = pg.LinearRegionItem([0.104, 0.112])
        for time_rgn in (self.base_time_rgn, self.test_time_rgn):
            self.plt1.addItem(time_rgn)
            time_rgn.sigRegionChangeFinished.connect(self.time_rgn_changed)

        self.plot_data = []

    def load_data(self, seq_dir):
        with pg.BusyCursor():
            self._load_data(seq_dir)

    def _load_data(self, seq_dir):
        self.clear_plot_data()

        man = getManager()
        model = man.dataModel

        # read all image data
        self.img_data = model.buildSequenceArray(
            seq_dir, 
            lambda dh: dh['Camera']['frames.ma'].read()['Time':0:200e-3].asarray(),
            join=True)
        
        first_subdir = seq_dir[seq_dir.ls()[0]]
        first_img = first_subdir['Camera']['frames.ma'].read()['Time':0:200e-3]
        self.img_data._info[2] = first_img._info[0]

        time_vals = self.img_data.xvals('Time')
        time_prof = self.img_data[:, 1].mean(axis=0).mean(axis=1).mean(axis=1)
        self.plot_data.append(self.plt1.plot(time_vals, time_prof))

        self.time_rgn_changed()

    def time_rgn_changed(self):
        base_start, base_stop = self.base_time_rgn.getRegion()
        test_start, test_stop = self.test_time_rgn.getRegion()
        base = self.img_data['Time':base_start:base_stop].asarray().mean(axis=2)
        test = self.img_data['Time':test_start:test_stop].asarray().mean(axis=2)

        base_diff = base[:,0] - base[:,1]
        test_diff = test[:,0] - test[:,1]

        self.diff_img = test_diff - base_diff
        self.imv1.setImage(ndimage.median_filter(self.diff_img, (3, 3, 3)))

    def clear_plot_data(self):
        for item in self.plot_data:
            self.plt1.removeItem(item)
        self.plot_data = []
