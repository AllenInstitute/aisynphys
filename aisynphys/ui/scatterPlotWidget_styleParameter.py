from pyqtgraph import QtWidgets, QtCore
from pyqtgraph import parametertree as ptree
import numpy as np
from collections import OrderedDict
from pyqtgraph import functions as fn

class StyleMapParameter(ptree.types.GroupParameter):
    sigStyleChanged = QtCore.Signal(object)
    
    def __init__(self):
        self.fields = {}
        ptree.types.GroupParameter.__init__(self, name='Style', addText='Add Style..', addList=[])
        self.sigTreeStateChanged.connect(self.styleChanged)

    def styleChanged(self):
        self.sigStyleChanged.emit(self)

    def fieldNames(self):
        return list(self.fields.keys())

    def setFields(self, fields):
        self.fields = OrderedDict(fields)
        names = self.fieldNames()
        self.setAddList(names)

    def addNew(self, name):
        fieldSpec = self.fields[name]
        
        mode = fieldSpec.get('mode', 'range')        
        if mode == 'range':
            item = RangeStyleItem(name, self.fields[name])
        elif mode == 'enum':
            item = EnumStyleItem(name, self.fields[name])
            
        self.addChild(item)
        return item

    def map(self, data):
        if isinstance(data, dict):
            data = np.array([tuple(data.values())], dtype=[(k, float) for k in data.keys()])

        symbols = np.full(len(data), 'o', dtype=str)
        # symbol_size = np.full(len(data), 5, dtype=int)
        # symbol_pen = dict(color=np.full(len(data), fn.mkBrush(1, 1, 1), dtype=object), width=np.full(len(data), 1., dtype=float))
        style = dict(pen=None, symbol=symbols)
        for item in self.children():
            if item.value is False:
                continue

            style = item.map(data)

        return style

class RangeStyleItem(ptree.types.SimpleParameter):
    mapType = 'range'

    def __init__(self, name, opts):
        self.fieldName = name
        units = opts.get('units', '')
        ptree.types.SimpleParameter.__init__(self, 
            name=name, autoIncrementName=True, removable=True, renamable=True, type='bool', value=True,
            children=[
                #dict(name="Field", type='list', value=name, values=fields),
                dict(name='Min', type='float', value=0.0, suffix=units, siPrefix=True),
                dict(name='Max', type='float', value=1.0, suffix=units, siPrefix=True),
                dict(name='Symbol', type='list', values=['o', 's', 't', 't1', '+', 'd'], value='o'),
                dict(name='Symbol size', type='int', value=10),
                dict(name='Symbol pen', type='group', expanded=False, children=[
                        dict(name='Color', type='color', value=fn.mkColor('w')),
                        dict(name='Width', type='float', value=1.0),
                    ])
            ])

    def map(self, data):
        vals = data[self.fieldName]
        if len(vals) == 0:
            return

        symbols = np.full(len(vals), 'o', dtype=str)
        symbol_size = np.full(len(data), 10, dtype=int)
        symbol_pen = np.full(len(data), fn.mkPen('w', width=1.0), dtype=object)
        mask = (vals >= self['Min']) & (vals < self['Max'])
        symbols[mask] = self['Symbol']
        symbol_size[mask] = self['Symbol size']
        symbol_pen[mask] = fn.mkPen(self['Symbol pen', 'Color'], width=self['Symbol pen', 'Width'])
        # symbol_pen['width'][mask] = self['Symbol pen', 'Width']
        
        style = dict(pen=None, symbol=symbols, symbolSize=symbol_size, symbolPen=symbol_pen)
        
        return style 

class EnumStyleItem(ptree.types.SimpleParameter):
    mapType = 'enum'

    def __init__(self, name, opts):
        self.fieldName = name
        vals = opts.get('values', [])
        if isinstance(vals, list):
            vals = OrderedDict([(v,str(v)) for v in vals])
        # childs = [{'name': v, 'type': 'group', 'children': [
        #             {'name': 'Symbol', 'type': 'list', 'values':['o', 's', 't', 't1', '+', 'd'], 'value': 'o'},
        #             {'name': 'Symbol size', 'type': 'int', 'value': 10},
        #             {'name': 'Symbol pen', 'type': 'group', 'expanded': False, 'children': [
        #                 {'name': 'Color', 'type': 'color', 'value': fn.mkColor('w')},
        #                 {'name': 'Width', 'type': 'float', 'value': 1.0},
        #             ]}
        # ]} for v in vals]
        
        childs = []
        for val,vname in vals.items():
            ch = ptree.Parameter.create(name=vname, type='bool', value=False, children=[
                dict(name='Symbol', type='list', values=['o', 's', 't', 't1', '+', 'd'], value='o'),
                dict(name='Symbol size', type='int', value=10),
                dict(name='Symbol pen', type='group', expanded=False, children=[
                        dict(name='Color', type='color', value=fn.mkColor('w')),
                        dict(name='Width', type='float', value=1.0),
                    ])
                ])
            ch.maskValue = val
            childs.append(ch)
           
        ptree.types.SimpleParameter.__init__(self, 
            name=name, autoIncrementName=True, removable=True, renamable=True, type='bool', value=True,
            children=childs)
                # dict(name='Values', type='group', children=childs),
                # dict(name='Symbol', type='list', values=['o', 's', 't', 't1', '+', 'd'], value='o'),
                # dict(name='Symbol size', type='int', value=10),
                # dict(name='Symbol pen', type='group', expanded=False, children=[
                #         dict(name='Color', type='color', value=fn.mkColor('w')),
                #         dict(name='Width', type='float', value=1.0),
                #     ])
            # ])
        
    def map(self, data):
        vals = data[self.fieldName]
        if len(vals) == 0:
            return

        symbols = np.full(len(vals), 'o', dtype=str)
        symbol_size = np.full(len(data), 10, dtype=int)
        symbol_pen = np.full(len(data), fn.mkPen('w', width=1.0), dtype=object)

        for child in self.children():
            if child.value() is False:
                continue
            mask = vals == child.maskValue
            symbols[mask] = child['Symbol']
            symbol_size[mask] = child['Symbol size']
            symbol_pen[mask] = fn.mkPen(child['Symbol pen', 'Color'], width=child['Symbol pen', 'Width'])

        return dict(pen=None, symbol=symbols, symbolSize=symbol_size, symbolPen=symbol_pen)
