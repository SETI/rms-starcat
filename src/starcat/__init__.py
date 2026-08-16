################################################################################
# starcat/__init__.py
################################################################################

from starcat.spice import SpiceStar, SpiceStarCatalog
from starcat.starcatalog import (
    SCLASS_TO_B_MINUS_V,
    SCLASS_TO_SURFACE_TEMP,
    Star,
    StarCatalog,
)
from starcat.ucac4 import UCAC4Star, UCAC4StarCatalog
from starcat.ybsc import YBSCStar, YBSCStarCatalog

# Version
try:
    from ._version import __version__
except ImportError:
    __version__ = 'Version unspecified'

__all__ = [
    'SCLASS_TO_B_MINUS_V',
    'SCLASS_TO_SURFACE_TEMP',
    'SpiceStar',
    'SpiceStarCatalog',
    'Star',
    'StarCatalog',
    'UCAC4Star',
    'UCAC4StarCatalog',
    'YBSCStar',
    'YBSCStarCatalog',
    '__version__'
]
