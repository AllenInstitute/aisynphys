
from __future__ import print_function, division

from collections import OrderedDict
from . import database as db
from . import constants


class CellClass(object):
    """Represents a class of cells as a list of selection criteria.
    
    Construct with an arbitrary set of keyword arguments, where each argument specifies
    a criteria for matching cells. Keyword argument names must be a column from the Cell
    or Morphology database tables.
    
    Examples
    --------
    
        pv_class = CellClass(cre_type='pvalb')
        inhibitory_class = CellClass(cre_type=('pvalb', 'sst', 'vip'))
        l23_pyr_class = CellClass(pyramidal=True, target_layer='2/3')
        l5_spiny_class = CellClass(dendrite_type='spiny', cortical_layer='5')
    """
    def __init__(self, display_names=None, **criteria):
        self.criteria = criteria
        self.display_names = display_names

    @property
    def name(self):
        """A short string representation of this cell class.
        
        This is the same as as_tuple, concatenated with spaces.
        """
        return ' '.join(self.as_tuple)

    @property
    def as_tuple(self):
        """A tuple representation of this cell class used for display purposes.
        
        Order of elements in the tuple is (target_layer, pyramidal, cre_type), but
        elements are only present if they were specified as criteria for the cell class.
        """
        if self.display_names is not None:
            return self.display_names
            
        name = []

        target_layer = self.criteria.get('target_layer')
        cortical_layer = self.criteria.get('cortical_layer')

        if target_layer is not None:
            name.append('L' + target_layer)
        elif cortical_layer is not None:
            name.append('L' + cortical_layer)

        if 'dendrite_type' in self.criteria:
            name.append('%s' % self.criteria['dendrite_type'])

        if 'pyramidal' in self.criteria:
            name.append('pyr' if self.criteria['pyramidal'] else 'nonpyr')

        cre_type = self.criteria.get('cre_type')
        if cre_type is not None:
            name.append(str(cre_type))
        
        return tuple(name)

    @property
    def is_excitatory(self):
        """True if this is an excitatory cell class, as determined either by
        cre type or morphology.
        """
        cre = self.criteria.get('cre_type')
        pyr = self.criteria.get('pyramidal')
        dendrite = self.criteria.get('dendrite_type')
        return cre == 'unknown' or cre in constants.EXCITATORY_CRE_TYPES or pyr is True or dendrite == 'spiny'

    def __contains__(self, cell):
        morpho = cell.morphology
        objs = [cell, morpho]
        for k, v in self.criteria.items():
            found_attr = False
            for obj in objs:
                if hasattr(obj, k):
                    found_attr = True
                    if isinstance(v, tuple):
                        if getattr(obj, k) not in v:
                            return False
                    else:
                        if getattr(obj, k) != v:
                            return False
                    break
            if not found_attr:
                raise Exception('Cannot use "%s" for cell typing; attribute not found on cell or cell.morphology' % k)
        return True

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, a):
        """Cell class is considered equal to its *name* to allow it to be indexed from a dict more
        easily::

            cc = CellClass(cre_type='sst', layer='6')
            cc.name => 'L6 sst'
            {cc: 1}['L6 sst'] => 1 
        """
        if isinstance(a, str):
            return a == self.name
        elif isinstance(a, CellClass):
            return a.name == self.name
        else:
            return object.__eq__(self)  # should raise NotImplemented

    def __repr__(self):
        return "<CellClass %s>" % self.name

    def __str__(self):
        return self.name
        
    def filter_query(self, query, cell_table):
        """Return a modified query (sqlalchemy) that filters results to include only those in
        this cell class.
        """
        morpho = db.aliased(db.Morphology)
        query = query.join(morpho, morpho.cell_id==cell_table.id)
        tables = [cell_table, morpho]
        for k, v in self.criteria.items():
            found_attr = False
            for table in tables:
                if hasattr(table, k):
                    found_attr = True
                    if isinstance(v, tuple):
                        query = query.filter(getattr(table, k).in_(v))
                    else:
                        query = query.filter(getattr(table, k)==v)
                    break
            if not found_attr:
                raise Exception('Cannot use "%s" for cell typing; attribute not found on cell or cell.morphology' % k)

        return query                
                
            


def classify_cells(cell_classes, cells=None, pairs=None, session=None):
    """Given cell class definitions and a list of cells, return a dict indicating which cells
    are members of each class.

    Parameters
    ----------
    cell_classes : dict
        List of CellClass instances
    cells : list | None
        List of Cell instances to be classified.
    pairs : list | None
        List of pairs from which cells will be collected. May not be used with *cells* or *session*
    session: Session | None
        If *cells* is not provided, then a database session may be given instead from which
        cells will be selected.
        
    Returns
    -------
    cell_groups : OrderedDict
        Dictionary mapping {cell_class: [list of cells]}
        
        
    Example
    -------
        
        pv_cell_class = CellClass(cre_type='pvalb', target_layer='2/3')
        sst_cell_class = CellClass(cre_type='sst', target_layer='2/3')
        
        cell_classes = [pv_cell_class, sst_cell_class]
        cells = session.Query(db.Cell).all()
        
        grouped_cells = classify_cells(cell_classes, cells=cells)
    
        pv_cells = grouped_cells[pv_cell_class]    
        sst_cells = grouped_cells[sst_cell_class]    

    """
    if pairs is not None:
        assert cells is None, "cells and pairs arguments are mutually exclusive"
        assert session is None, "session and pairs arguments are mutually exclusive"
        cells = set([p.pre_cell for p in pairs] + [p.post_cell for p in pairs])
    if cells is None:
        cells = session.query(db.Cell, db.Cell.cre_type, db.Cell.target_layer, db.Morphology.pyramidal, db.Morphology.cortical_layer, db.Morphology.dendrite_type).join(db.Morphology)
    cell_groups = OrderedDict([(cell_class, set()) for cell_class in cell_classes])
    for cell in cells:
        for cell_class in cell_classes:
            if cell in cell_class:
                cell_groups[cell_class].add(cell)
    return cell_groups


def classify_pairs(pairs, cell_groups):
    """Given a list of cell pairs and a dict that groups cells together by class (ie the output of classify_cells),
    return a dict that groups pairs into (pre, post) cell type buckets.
    
    Parameters
    ----------
    pairs : list of Pair instances
        The Pair instances (probably returned from a database query) to be grouped
    cell_groups : dict
        Specifies the cell classes and the cells that belong to each class. The format is
        the same as the output of classify_cells().
        
    Returns
    -------
    pair_groups : OrderedDict
        Maps {(pre_class, post_class): [list of pairs]}
    """
    results = OrderedDict()
    for pre_class, pre_group in cell_groups.items():
        for post_class, post_group in cell_groups.items():
            post_group = cell_groups[post_class]
            class_pairs = [p for p in pairs if p.pre_cell in pre_group and p.post_cell in post_group]
            results[(pre_class, post_class)] = class_pairs
    
    return results
