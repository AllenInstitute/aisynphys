import sys

def interactive_mode():
    """Return a string describing the preferred mode of user interaction.
    
    Can be one of: 'qt', 'ipynb', 'tty', or 'file'.
    """
    if 'pyqtgraph' in sys.modules and sys.modules['pyqtgraph'].QtWidgets.QApplication.instance() is not None:
        return 'qt'

    if 'IPython' in sys.modules:
        kern = sys.modules['IPython'].get_ipython()
        if 'ZMQ' in str(kern):
            return 'ipynb'
        else:
            return 'tty'

    if sys.stdout.isatty():
        return 'tty'

    return 'file'
