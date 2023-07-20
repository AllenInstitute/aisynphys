from .genotypes import FLUOROPHORES

INHIBITORY_CRE_TYPES = ['sst', 'pvalb', 'vip', 'ndnf', 'chat', 'htr3a', 'nos1', 'chrna2', 'mDlx', '3xcorehI56i', 'tac1']
EXCITATORY_CRE_TYPES = ['tlx3', 'sim1', 'fam84b', 'rorb', 'ntsr1', 'rbp4', 'ctgf', 'glt25d2', 'slc17a8', 'cux2', 'nr5a1', 'eHGT_078m', '3xHGT_073m(core)']
ALL_CRE_TYPES = INHIBITORY_CRE_TYPES + EXCITATORY_CRE_TYPES + ['unknown']

DRIVER_LAYERS = {
    'cux2':    ['2', '2/3', '3', '4'],
    'rorb':    ['4'],
    'nr5a1':   ['4'],
    'tlx3':    ['5'],
    'sim1':    ['5'],
    'rbp4':    ['5'],
    'slc17a8': ['5', '6'],
    'ntsr1':   ['6'],
    'ctgf':    ['6'],
}

# note: the lower-case af488 and cascade_blue are for backward compatibility; these
# may be removed at some point
ALL_LABELS = ['biocytin', 'af488', 'af594', 'cascade_blue'] + list(FLUOROPHORES.keys())

LAYERS = ['1', '2', '2/3', '3', '4', '5', '5a', '5b', '6']

ACSF_RECIPES = ["Standard", "ACSF III (no blockers)", "2mM Ca & Mg", "1.3mM Ca & 1mM Mg", "4mM Ca", "4mM Ca + TTx,4AP", "mixed"]
INTERNAL_RECIPES = ["Standard K-Gluc", "K-Gluc -EGTA", "K-Gluc 1uM EGTA", 'PatchSeq', 'Cesium', 'mixed']
INTERNAL_DYES = ['Cascade Blue', 'AF488', 'AF594', 'no dye']


INJECTIONS = {
    'RO Fam84b-FlpO': 'Fam84b-FlpO',
    'RO tdTomato+GFP': 'pAAV-Ef1a-fDIO-EGFP;pAAV-Ef1a-cDIO-dTomato',
    'pan-GABA': 'rAAV-Dlx2.0-SYFP2',
    'CN1390': 'rAAV-3xhI56i(core)-minBglobin-SYFP2-WPRE3-BGHpA',
    'CN1461': 'rAAV-eHGT_073m-minBglobin-SYFP2-WPRE3-BGHpA',
    'CN1466': 'rAAV-eHGT_078m-minBglobin-SYFP2-WPRE3-BGHpA',
    'CN1809': 'rAAV-TRE-tdTomato-WPRE-HGHpA',
    'CN1810': 'rAAV-TRE-SYFP2-WPRE-HGHpA',
    'CN1821': 'rAAV-hSyn1-tTA-WPRE-HGHpA',
    'CN1827': 'rAAV-3xhI56icore-minBG-tdTomato-WPRE3-BGHpA',
    'CN1849': 'rAAV-3xHGT_073m(core)-minBG-SYFP2-WPRE3-BGHpA',
    'CN1913': 'rAAV-eHGT_078m-minBG-FlpO-WPRE-HGHpA',
    'CN1914': 'rAAV-eHGT_078m-minBG-tTA-WPRE-HGHpA',
    'CN1915': 'rAAV-eHGT_073h-minBG-tTA-WPRE-HGHpA',
    'CN1955': 'rAAV-hsA2-3x(eHGT_078m core)-minRho-SYFP2-WPRE3-BGHpA',
    'CN1987': 'rAAV-3xeHGT_073m(core)-minCMV*-SYFP2-WPRE3-BGHpA',
    'CN1988': 'rAAV-EF1a-fDIO-EGFP-WPRE-HGHpA (Miranda/Tanya)',
    'CN1995': 'rAAV-TREtight-tdTomato-WPRE-HGHpA',
    'ChrimsonR-Kv-EYFP':'AAV-Syn-FLEX-ChrimsonR-eYFP-Kv2.1',
}
