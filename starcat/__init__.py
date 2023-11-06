################################################################################
# starcat/__init__.py
################################################################################

from starcat.spice import *
from starcat.ucac4 import *
from starcat.ybsc import *

try:
    from _version import __version__
except ImportError as err:
    __version__ = 'Version unspecified'
