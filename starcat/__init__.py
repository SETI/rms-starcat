################################################################################
# starcat/__init__.py
################################################################################

from starcat.starcatalog import Star, StarCatalog  # noqa: F401
from starcat.spice import SpiceStar, SpiceStarCatalog  # noqa: F401
from starcat.ucac4 import UCAC4Star, UCAC4StarCatalog  # noqa: F401
from starcat.ybsc import YBSCStar, YBSCStarCatalog  # noqa: F401

# Version
try:
    from ._version import __version__
except ImportError:
    __version__ = 'Version unspecified'

__all__ = [
    'Star',
    'StarCatalog',
    'SpiceStar',
    'SpiceStarCatalog',
    'UCAC4Star',
    'UCAC4StarCatalog',
    'YBSCStar',
    'YBSCStarCatalog',
    '__version__'
]
