import sys

def interactive_mode():
    """Return a string describing the preferred mode of user interaction.
    
    Can be one of: 'qt', 'ipynb', 'tty', or 'file'.
    """
    try:
        if 'pyqtgraph' in sys.modules:
            pg = sys.modules['pyqtgraph']
            if hasattr(pg, 'QtWidgets') and pg.QtWidgets.QApplication.instance() is not None:
                return 'qt'
            elif hasattr(pg, 'Qt') and pg.Qt.QtWidgets.QApplication.instance() is not None:
                return 'qt'
    except Exception:
        pass

    if 'IPython' in sys.modules:
        kern = sys.modules['IPython'].get_ipython()
        if 'ZMQ' in str(kern):
            return 'ipynb'
        else:
            return 'tty'

    if sys.stdout.isatty():
        return 'tty'

    return 'file'
