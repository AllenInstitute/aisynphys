from __future__ import print_function
import sys, itertools

# ----------------- Genotype Handling ----------------------------
# The purpose of the code below is to be able to ask, given a mouse genotype string, 
# what cell types are marked, and what reporter color are they marked with.

# There are two approaches:
# 1. Write down a list of all possible genotype strings and the associated 
#    driver / reporter information. This is simple and transparent, but requires
#    constant maintenance as we add more genotypes.
# 2. Write down the complete list of drivers and reporters, and use some simple
#    genetic modeling to predict what phenotype will be created by a particular
#    genotype.
# The Genotype class below implements the second approach.


# Background research on understanding genotype symbols:
# http://help.brain-map.org/download/attachments/2818171/Connectivity_Resources.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4365051/pdf/nihms665520.pdf

# RCL  : Rosa_CAG-LSL      (cre dependent)
# RCFL : Rosa_CAG-FSF-LSL  (cre+flp dependent)
# RCRL : Rosa_CAG-RSR-LSL  (cre+dre dependent)
# TIT : TIGRE-insulators-TRE_promoter   (tTa|rtTA dependent)
# TITL : TIGRE-insulators-TRE_promoter-LSL   (cre+(tTA|rtTA) dependent)
# TIT2L : TIGRE-insulators-TRE_promoter-LSL   (cre+(tTA|rtTA) dependent)
# IRES :
# T2A :
# LSL : lox-stop-lox (cre dependent)
# FSF : frt-stop-frt (flp dependent)
# RSR : rox-stop-rox (dre dependent, but also slightly activated by cre)
# dgCre : functions like Cre with reduced activity; increased by administration of TMP
# TRE : tetracycline response element - series of TetO sequences that
#       require tTA or rtTA to increase downstream expression
#         - dox binds tTa, which prevents further binding to TetO
#         - dox binds rtTa, which allows further binding to TetO
#       this system is also used for amplification.
# tTA : tet-off  (activates TRE only in absence of doxycycline)
# rtTA: tet-on   (activates TRE only in presence of doxycycline)
# tTA2 : same behavior as tTA

# hI56i (or hI56(1)) : human, intergenic region between dlx5 and dlx6, enhancer i  (i is pan-GABA; ii is not)
#                      3xcore : minimal (bashed) core region of enhancer, repeated 3x

EXPRESSION_FACTORS = ['cre', 'flp', 'dre', 'tTA']
DRUGS = ['dox']

   
DRIVER_LINES = {                  # dependencies     products
    'Nr5a1-Cre':                  [(['nr5a1'],        ['cre'])],
    'Rorb-T2A-tTA2':              [(['rorb'],         ['tTA'])],
    'Tlx3-Cre_PL56':              [(['tlx3'],         ['cre'])],
    'Sim1-Cre_KJ18':              [(['sim1'],         ['cre'])],
    'Ntsr1-Cre_GN220':            [(['ntsr1'],        ['cre'])],
    'Chat-IRES-Cre-neo':          [(['chat'],         ['cre'])],
    'Rbp4-Cre_KL100':             [(['rbp4'],         ['cre'])],
    'Pvalb-IRES-Cre':             [(['pvalb'],        ['cre'])],
    'Sst-IRES-Cre':               [(['sst'],          ['cre'])],
    'Vip-IRES-Cre':               [(['vip'],          ['cre'])],
    'Pvalb-T2A-FlpO':             [(['pvalb'],        ['flp'])],
    'Sst-IRES-FlpO':              [(['sst'],          ['flp'])],
    'Vip-IRES-FlpO':              [(['vip'],          ['flp'])],
    'Slc17a8-IRES2-Cre':          [(['slc17a8'],      ['cre'])],
    'Pvalb-2A-FlpO':              [(['pvalb'],        ['flp'])],
    'Cux2-CreERT2':               [(['cux2'],         ['cre'])],
    'Chrna2-Cre_OE25':            [(['chrna2'],       ['cre'])],
    'Penk-IRES2-Cre-neo':         [(['penk'],         ['cre'])],
    'Slc17a6-IRES2-FlpO':         [(['slc17a6'],      ['flp'])],
    'Slc17a8-iCre':               [(['slc17a8'],      ['cre'])],
    'Ctgf-T2A-dgCre':             [(['ctgf'],         ['cre'])],
    'Ndnf-IRES2-dgCre':           [(['ndnf'],         ['cre'])],
    'Slc32a1-IRES2-FlpO':         [(['slc32a1'],      ['flp'])],
    'Fam84b-FlpO':                [(['fam84b'],       ['flp'])],
    'rAAV-mDlx-GFP':              [(['mDlx'],         ['GFP'])],    
    'rAAV-Dlx2.0-SYFP2':          [(['3xcorehI56i'],  ['YFP'])],  # pan-GABA
    'rAAV-eHGT_078m-minBglobin-SYFP2-WPRE3-BGHpA': [
                                   (['eHGT_078m'], ['YFP'])],  # pan-Glu
    'rAAV-3xhI56icore-minBG-tdTomato-WPRE3-BGHpA': [
                                   (['3xcorehI56i'],  ['tdTomato'])],  # pan-GABA
    'rAAV-eHGT_073m-minBglobin-SYFP2-WPRE3-BGHpA': [
                                   (['eHGT_073m'], ['YFP'])],  # pan-Glu
    'rAAV-3xHGT_073m(core)-minBG-SYFP2-WPRE3-BGHpA': [
                                   (['3xHGT_073m(core)'], ['YFP'])],  # pan-Glu
    'rAAV-EF1a-fDIO-EGFP-WPRE-HGHpA': [
                                   (['EF1a'], ['EGFP'])],
    'rAAV-eHGT_078m-minBG-FlpO-WPRE-HGHpA': [
                                   (['eHGT_078m'], ['flp'])],  # pan-Glu
}


REPORTER_LINES = {                # dependencies             products
    'Ai2':                        [([],                       ['EYFP'])],
    'Ai3':                        [([],                       ['EYFP'])],
    'Ai6':                        [([],                       ['ZsGreen'])],
    'Ai9':                        [([],                       ['tdTomato'])],
    'Ai14':                       [([],                       ['tdTomato'])],
    'Ai14(RCL-tdT)':              [(['cre'],                  ['tdTomato'])],
    'Ai27':                       [([],                       ['hChR2(H134R)', 'tdTomato'])],
    'Ai31':                       [([],                       ['Syp', 'Emerald'])],
    'Ai32':                       [([],                       ['ChR2(H134R)', 'EYFP'])],
    'Ai34':                       [([],                       ['Syp', 'tdTomato'])],
    'Ai35':                       [([],                       ['Arch', 'EGFP', 'ER2'])],
    'Ai57(RCFL-Jaws)':            [(['cre', 'flp'],           ['Jaws', 'GFP', 'ER2'])],
    'Ai62(TITL-tdT)':             [(['cre', 'tTA', '~dox'],   ['tdTomato']),
                                   (['cre', 'rtTA', 'dox'],   ['tdTomato'])],
    'Ai63(TIT-tdT)':              [(['tTA', '~dox'],          ['tdTomato']),
                                   (['rtTA', 'dox'],          ['tdTomato'])],
    'Ai65(RCFL-tdT)':             [(['cre', 'flp'],           ['tdTomato'])],
    'Ai65F':                      [(['flp'],                  ['tdTomato'])],
    'Ai66(RCRL-tdT)':             [(['cre', 'dre'],           ['tdTomato'])],
    'Ai72(RCL-VSFPB)':            [(['cre'],                  ['VSFP', 'Butterfly 1.2'])],
    'Ai78(TITL-VSFPB)':           [(['cre', 'tTA', '~dox'],   ['VSFP', 'Butterfly 1.2']),
                                   (['cre', 'rtTA', 'dox'],   ['VSFP', 'Butterfly 1.2'])],
    'Ai79(TITL-Jaws)':            [(['cre', 'tTA', '~dox'],   ['Jaws', 'GFP', 'ER2']),
                                   (['cre', 'rtTA', 'dox'],   ['Jaws', 'GFP', 'ER2'])],
    'Ai82(TITL-GFP)':             [(['cre', 'tTA', '~dox'],   ['EGFP']),
                                   (['cre', 'rtTA', 'dox'],   ['EGFP'])],
    'Ai85(TITL-iGluSnFR)':        [(['cre', 'tTA', '~dox'],   ['iGluSnFR']),
                                   (['cre', 'rtTA', 'dox'],   ['iGluSnFR'])],
    'Ai87(RCL-iGluSnFR)':         [(['cre'],                  ['iGluSnFR'])],
    'Ai92(TITL-YCX2.60)':         [(['cre', 'tTA', '~dox'],   ['YCX2.60']),
                                   (['cre', 'rtTA', 'dox'],   ['YCX2.60'])],
    'Ai93(TITL-GCaMP6f)':         [(['cre', 'tTA', '~dox'],   ['GCaMP6f']),
                                   (['cre', 'rtTA', 'dox'],   ['GCaMP6f'])],
    'Ai94(TITL-GCaMP6s)':         [(['cre', 'tTA', '~dox'],   ['GCaMP6s']),
                                   (['cre', 'rtTA', 'dox'],   ['GCaMP6s'])],
    'Ai95(RCL-GCaMP6f)':          [(['cre'],                  ['GCaMP6f'])],
    'Ai96(RCL-GCaMP6s)':          [(['cre'],                  ['GCaMP6s'])],
    'Ai139(TIT2L-GFP-ICL-TPT)-D': [(['cre'],                  ['EGFP', 'tdTomato'])],
    'Ai139(TIT2L-GFP-ICL-TPT)':   [(['cre'],                  ['EGFP', 'tdTomato'])],
    'Ai140(TIT2L-GFP-ICL-tTA2)':  [(['cre'],                  ['EGFP', 'tTA'])],
    'Ai167(TIT2L-ChrimsonR-tdT-ICL-tTA2)': [
                                   (['cre'],                  ['ChrimsonR', 'tdTomato'])],
    'Snap25-LSL-F2A-GFP':         [(['cre'],                  ['EGFP'])],
    'Ai193-hyg-440167':           [(['cre'],                  ['EGFP']),
                                   (['flp'],                  ['tdTomato'])],
    'Ai193(TICL-EGFP-ICF-tdT)-hyg': [
                                   (['cre'],                  ['EGFP']),
                                   (['flp'],                  ['tdTomato'])],
    'pAAV-Ef1a-fDIO-EGFP':        [(['flp'],                  ['EGFP'])],
    'pAAV-Ef1a-cDIO-dTomato':     [(['cre'],                  ['tdTomato'])],
}


FLUOROPHORES = {
    'tdTomato': 'red',
    'GFP': 'green',
    'EGFP': 'green',
    'AF488': 'green',
    'Cascade Blue': 'blue',
    'YFP': 'yellow',
    'EYFP': 'yellow',
    'ZsGreen': 'green',
    'GCamp6f': 'green',
    'Emerald': 'green',
}


ALL_TG_LINES = DRIVER_LINES.copy()
ALL_TG_LINES.update(REPORTER_LINES)

ALL_COLORS = set(FLUOROPHORES.values())

ALL_DRIVERS = set()
for rules in DRIVER_LINES.values():
    for deps, prods in rules:
        ALL_DRIVERS |= set([d for d in deps if d[0] != '~'])
ALL_DRIVERS = ALL_DRIVERS - set(EXPRESSION_FACTORS) - set(DRUGS)

ALL_REPORTERS = set(FLUOROPHORES.keys())
for rules in REPORTER_LINES.values():
    for deps, prods in rules:
        ALL_REPORTERS |= set(prods)
ALL_REPORTERS = ALL_REPORTERS - set(EXPRESSION_FACTORS) - set(DRUGS)



class Genotype(object):
    """
    Class used to provide phenotype information from a genotype string.

    Parameters
    ----------
    gtype : str
        A genotype string like
        'Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt'

    Notes
    -----

    In the example given above, the genotype specifies two driver lines (Tlx3-Cre_PL65
    and Sst-IRES-FlpO) and two reporter lines (Ai65F and Ai140). This results in
    tdTomato expression in sst cells and GFP expression in tlx3 cells.

    Example usage::

        gtype_str = 'Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt'
        gtype = Genotype(gtype_str)

        # Parsed genotype parts
        gtype.driver_lines                       # =>  ['Sst-IRES-FlpO', 'Tlx3-Cre_PL56']
        gtype.reporter_lines                     # =>  ['Ai140(TIT2L-GFP-ICL-tTA2)', 'Ai65F']

        # Reduced representation of genotype parts
        gtype.all_drivers                        # =>  set(['sst', 'tlx3'])
        gtype.all_reporters                      # =>  set(['EGFP', 'tdTomato'])
        gtype.all_colors                         # =>  set(['green', 'red'])

        # What reporters are expressed if the cell is sst positive?
        gtype.expressed_reporters(['tlx3'])      # =>  set(['EGFP'])

        # What colors are expressed if the cell is sst positive?
        gtype.expressed_colors(['sst'])          # =>  set(['red'])

        # Given colors seen in a cell, can we predict which drivers were active?
        colors = {'red': True, 'green': False}
        gtype.predict_driver_expression(colors)  # =>  {'sst': True, 'tlx3': False}

    """
    def __init__(self, gtype):
        self.gtype = gtype
        self._parse()

    def __repr__(self):
        return "<Genotype %s>" % self.gtype
    
    def expressed_reporters(self, drivers):
        """Return a set of reporters in this genotype (such as GFP, tdTomato, etc.)
        that would be expressed given a list of active drivers (such as pvalb, tlx3, etc.)

        Parameters
        ----------
        drivers : set
            The set of active drivers

        Notes
        -----

        Example::

            # A simple case quad-transgenic genotype where tlx3->EGFP, pvalb->tdTomato
            gt = Genotype('Tlx3-Cre_PL56/wt;Pvalb-2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt')

            gt.expressed_reporters(['pvalb'])
            # returns: set(['tdTomato'])
        """
        prods = self.model.forward_model(drivers)
        reporters = set([p for p in prods if p in ALL_REPORTERS])
        return reporters

    def expressed_colors(self, drivers):
        """Return a set of fluorescent emission colors generated by reporters in this genotype
        that would be expressed given a list of active drivers (such as pvalb, tlx3, etc.)

        Parameters
        ----------
        drivers : set
            The set of active drivers

        Notes
        -----

        Example::

            # A simple case quad-transgenic genotype where tlx3->EGFP, pvalb->tdTomato
            gt = Genotype('Tlx3-Cre_PL56/wt;Pvalb-2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt')

            gt.expressed_colors(['pvalb'])
            # returns: set(['red'])
        """
        prods = self.model.forward_model(drivers)
        colors = set([p for p in prods if p in ALL_COLORS])
        return colors

    def predict_driver_expression(self, colors, starting_factors=()):
        """Given information about fluorescent colors expressed in a cell,
        return predictions about whether each driver could have been active.

        Parameters
        ----------
        colors : dict
            Describes whether each color observed in a cell was expressed (True),
            not expressed (False), or ambiguous (None).
        starting_factors : list
            Optional list of factors that are present (for example: dox)

        Returns
        -------
        driver_expression : dict
            Dict indicating whether each driver in the genotype must be expressed
            (True), must not be expressed (False), or cannot be determined (None).

        Notes
        -----

        Example::

            # A simple case quad-transgenic genotype where tlx3->EGFP, pvalb->tdTomato
            gt = Genotype('Tlx3-Cre_PL56/wt;Pvalb-2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt')

            gt.predict_driver_expression({'red': False, 'green': True})
            # returns: {'pvalb': False, 'tlx3': True}

            gt.predict_driver_expression({'red': True})
            # returns: {'pvalb': True, 'tlx3': None}
        """
        return self.model.reverse_model(unknown_factors=self.all_drivers, products=colors, starting_factors=starting_factors)

    def test_driver_combinations(self, colors, starting_factors=()):
        """Given information about fluorescent colors expressed in a cell,
        return predictions about whether each combination of drivers could
        possibly have generated the specified color measurements.

        Parameters
        ----------
        colors : dict
            Describes whether each color observed in a cell was expressed (True),
            not expressed (False), or ambiguous (None).
        starting_factors : list
            Optional list of factors that are present (for example: dox)
        
        Returns
        -------
        driver_expression : dict
            Dict indicating whether each driver combination in the genotype produces colors
            that are consistent (True) or inconsistent (False) with the observed colors.

        Notes
        -----

        The return values are only False where the predicted color expression is *inconsistent* with
        observed colors; in cases where there is not enough observed information, the return value
        will be True.

        Example::

            # A messy quad-transgenic genotype where pv->EGFP+tdTomato, rorb->tdTomato
            gt = Genotype('Pvalb-IRES-Cre/wt;Rorb-T2A-tTA2/wt;Ai63(TIT-tdT)/Ai140(TIT2L-GFP-ICL-tTA2)')  

            gt.test_driver_combinations({'red': True, 'green': False})
            # returns: {('pvalb',): False, ('rorb',): True, (): False, ('pvalb', 'rorb'): False}
            # this says that of the 4 possible combinations of drivers, only rorb+/pvalb- could 
            # generate red without green

            gt.test_driver_combinations({'red': True, 'green': True})
            # returns: {('pvalb',): True, ('rorb',): False, (): False, ('pvalb', 'rorb'): True}
            # this says that of the color expression is ambiguous--it could be the result of 
            # either rorb+/pvalb+ or rorb-/pvalb+


        """
        return self.model.test_factor_combinations(unknown_factors=self.all_drivers, products=colors, starting_factors=starting_factors)

    def _parse(self):
        """Extract driver/reporter lines from genotype string, generate a GeneticModel
        """
        ignore = ['wt', 'PhiC31-neo']
        
        parts = set()
        for part in self.gtype.split(';'):
            for subpart in part.split('/'):
                if subpart in ignore:
                    continue
                parts.add(subpart)
                
        self.driver_lines = [p for p in parts if p in DRIVER_LINES]
        self.reporter_lines = [p for p in parts if p in REPORTER_LINES]
        extra = parts - set(self.driver_lines + self.reporter_lines)
        if len(extra) > 0:
            raise Exception("Unknown genotype part(s): %s" % str(extra))

        # convert drivers / reporters into a GeneticModel
        self.model = GeneticModel()
        self.all_drivers = set()
        self.all_reporters = set()
        for line in self.driver_lines + self.reporter_lines:
            rules = ALL_TG_LINES[line]
            for deps, prods in rules:
                self.model.add_rule(deps, prods)
                self.all_drivers |= set([d.lstrip('~') for d in deps if d in ALL_DRIVERS])
                reporters = set([p for p in prods if p in ALL_REPORTERS])
                self.all_reporters |= reporters
                for reporter in reporters:
                    color = FLUOROPHORES.get(reporter, None)
                    if color is not None:
                        self.model.add_rule([reporter], [color])

        self.all_colors = set([FLUOROPHORES[r] for r in self.all_reporters if r in FLUOROPHORES])


class GeneticModel:
    """Genetic modeling engine.

    See add_rule(), forward_model(), and reverse_model().
    """

    def __init__(self, ruleset=None):
        # ruleset is a list of (inputs, outputs) tuples saying "if we have everything in inputs, then we generate everything in outputs"
        # e.g.:  inputs=('tlx3',)  outputs=('cre', 'tTA')
        #        inputs=('tTA',)   outputs=('tdTomato')
        self.ruleset = []
        self.all_products = set()
        if ruleset is not None:
            for dependencies, products in ruleset:
                self.add_rule(dependencies, products)

    def add_rule(self, dependencies, products):
        """Add a rule to the model that describes the transformation of a list of inputs to a list of outputs.

        For example, in a tlx3-cre animal the rule we want looks like "if cell expresses tlx3 then cell produces cre".
        That relationship would be written as::

            model.add_rule(dependencies=['tlx3'], products=['cre'])

        A rule can have any number of dependencies and products::

            # If cre AND flp are present, then produce EGFP AND ChR2
            model.add_rule(['cre', 'flp'], ['EGFP', 'ChR2'])

        Rules can also have negative dependencies::

            # If cell expresses tlx3 AND NOT sim1, then produce cre
            model.add_rule(['tlx3', '~sim1'], ['cre'])
        """
        assert isinstance(dependencies, (list, tuple, set)), "dependencies must be a list of strings"
        for dep in dependencies:
            assert isinstance(dep, str), "dependencies must be a list of strings"
        assert isinstance(products, (list, tuple, set)), "products must be a list of strings"
        for prod in products:
            assert isinstance(prod, str), "products must be a list of strings"
        pos_deps = set([d for d in dependencies if not d.startswith('~')])
        neg_deps = set([d[1:] for d in dependencies if d.startswith('~')])
        products = set(products)
        self.ruleset.append(((pos_deps, neg_deps), products))
        self.all_products |= products

    def forward_model(self, starting_factors):
        """Given a list of starting factors, predict the set of ending factors given the ruleset
        in this model.

        Example::

            # Imagine a simple genotype: tlx3->cre, cre->tdTomato
            model = GeneticModel()
            model.add_rule(['tlx3'], ['cre'])
            model.add_rule(['cre'], ['tdTomato'])

            # If a cell expresses tlx3, then it should also produce cre and tdTomato:
            model.forward_model(['tlx3'])   => set(['tlx3', 'cre', 'tdTomato'])

            # If a cell does not express tlx3, then nothing is produced
            model.forward_model([])   => set([])

        """
        factors = set(starting_factors)
        new_expression = False
        for dependencies, products in self.ruleset:
            pos_deps, neg_deps = dependencies

            # are all positive dependencies present, and no negative dependencies present?
            dependencies_met = pos_deps.issubset(factors) and len(factors & neg_deps) == 0

            # If dependencies are met, will any new products be generated?
            if dependencies_met and not products.issubset(factors):
                # Something new was expressed; add it to the list and run again
                factors |= products
                new_expression = True
        
        if new_expression:
            return self.forward_model(factors)
        else:
            return factors

    def reverse_model(self, unknown_factors, products, starting_factors=()):
        """Given information about products expressed in a cell,
        return predictions about whether each unknown starting factor could have been active.

        Parameters
        ----------
        unknown_factors : list
            List of factors may or may not have been present in the cell; the return value
            indicates for each of these products whether it was expressed.
        products : dict
            Describes whether each product measured in a cell was expressed (True),
            not expressed (False), or ambiguous (None).
        starting_factors : list
            List of factors that are known to exist

        Returns
        -------
        factor_expression : dict
            Dict indicating whether each unknown factor must be expressed
            (True), must not be expressed (False), or cannot be determined (None).

        Notes
        -----

        Example::

            # Imagine a quad-transgenic with the following relationships:
            #    tlx3 -> cre
            #    pvalb -> flp
            #    cre -> EGFP
            #    flp -> tdTomato
            model = GeneticModel()
            model.add_rule(['tlx3'], ['cre'])
            model.add_rule(['pvalb'], ['flp'])
            model.add_rule(['cre'], ['EGFP'])
            model.add_rule(['flp'], ['tdTomato'])

            # We observe a green fluorescent cell and want to know which factors
            # must have been expressed:
            unknown_factors=['tlx3', 'pvalb']
            products={'EGFP': True, 'tdTomato': False}
            model.reverse_model(unknown_factors, products)
            # returns: {'tlx3': True, 'pvalb': False}

        """
        factor_combos = self.test_factor_combinations(unknown_factors, products, starting_factors)
        true_combos = [factors for factors,prediction in factor_combos.items() if prediction is True]
        factor_expression = {}
        for factor in unknown_factors:
            n_match_with_factor = len([d for d in true_combos if factor in d])
            n_match_without_factor = len([d for d in true_combos if factor not in d])

            if n_match_with_factor > 0 and n_match_without_factor == 0:
                # If every possible factor combination that is consistent with observed products
                # contains this factor, then we say it is definitely expressed.
                factor_expression[factor] = True
            elif n_match_with_factor == 0:
                # If no possible combinations that contain this factor are consistent with observed
                # products, then we say the factor is definitely not expressed.
                # TODO: this can give incorrect results in cases where the factor is expressed, but silenced.
                factor_expression[factor] = False
            else:
                # Otherwise, we can't say one way or another whether this factor is expressed.
                factor_expression[factor] = None

        return factor_expression

    def test_factor_combinations(self, unknown_factors, products, starting_factors=()):
        """Given information about products expressed in a cell,
        return predictions about whether each combination of drivers could
        possibly have generated the products.

        Parameters
        ----------
        unknown_factors : list
            List of _possible_ starting factors to test. All combinations of these
            factors will be forward-modeled and compared against the observed products.
        products : dict
            Describes whether each product measured in a cell was expressed (True),
            not expressed (False), or ambiguous (None).
        starting_factors : list
            List of known starting factors to include during testing.
        
        Returns
        -------
        factor_expression : dict
            Dict indicating whether each starting factor combination produces products
            that are consistent (True) or inconsistent (False) with the observed colors.

        Notes
        -----

        The return values are only False where the predicted product expression is *inconsistent* with
        observed products; in cases where there is not enough observed information, the return value
        will be True.

            # Imagine a quad-transgenic with the following relationships:
            #    tlx3 -> cre
            #    pvalb -> flp
            #    cre -> EGFP
            #    flp -> tdTomato
            model = GeneticModel()
            model.add_rule(['tlx3'], ['cre'])
            model.add_rule(['pvalb'], ['flp'])
            model.add_rule(['cre'], ['EGFP'])
            model.add_rule(['flp'], ['tdTomato'])

            # We observe a green fluorescent cell and want to know what combination(s) of
            # starting factors (tlx3/pvalb) could have produced this color:
            unknown_factors=['tlx3', 'pvalb']
            products={'EGFP': True, 'tdTomato': False}
            model.test_factor_combinations(unknown_factors, products)
            # returns: {
            #    (): False,
            #    ('tlx3',): True,
            #    ('pvalb',): False,
            #    ('tlx3', 'pvalb'): False,
            # }

            Now try again, but this time we only observe EGFP and neglect to check for tdTomato:
            products={'EGFP': True}
            model.test_factor_combinations(unknown_factors, products)
            # returns: {
            #    (): False,
            #    ('tlx3',): True,
            #    ('pvalb',): False,
            #    ('tlx3', 'pvalb'): True,
            # }

        """
        starting_factors = set(starting_factors)
        predictions = {}
        # iterate over all combinations of unknown_factors
        for factors in self._factor_combinations(unknown_factors):
            # initial assumption is that this combination _will_ generate products that are
            # consistent with the observed products
            factor_combo_possible = True

            # model the products geenrated by this combination of factors
            predicted_products = self.forward_model(set(factors) | starting_factors)

            # check over all products for any mismatch
            for product in self.all_products:
                product_observed = products.get(product, None)
                if product_observed is None:
                    continue
                product_predicted = product in predicted_products
                if product_predicted != product_observed:
                    # prediction and observation do not match; mark this combination as impossible
                    factor_combo_possible = False
                    break
            predictions[tuple(sorted(factors))] = factor_combo_possible
        return predictions

    def _factor_combinations(self, factors):
        """Return a list of all possible combinations of the given factors"""
        factor_combos = []
        for i in range(len(factors) + 1):
            factor_combos.extend(list(itertools.combinations(factors, i)))
        return factor_combos

    def _simulate_factor_combos(self, starting_factors):
        """Given a set of starting factors that may or may not be present, forward-model 
        the products generated by all possible combinations of these factors.

        Returns
        -------
        driver_reporter_map : dict
            Maps {(drivers,): [reporters]} to describe all of the reporters that would be expressed
            by each possible combination of drivers.
        """
        product_map = {}

        # Iterate over all combinations of starting factors
        # predict which products would be expressed in each case
        factor_combos = self._factor_combinations(starting_factors)
        for factors in factor_combos:
            product_map[tuple(sorted(factors))] = self.forward_model(factors)

        return product_map
