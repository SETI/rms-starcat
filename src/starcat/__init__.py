################################################################################
# starcat/__init__.py
################################################################################

"""Read and search star catalogs.

This package provides a common interface to several star catalogs. Every catalog is
represented by a subclass of :class:`~starcat.StarCatalog` that yields stars as a
subclass of :class:`~starcat.Star`:

- :class:`~starcat.SpiceStarCatalog` for NAIF SPICE type 1 star catalogs, such as
  Hipparcos, PPM, and Tycho-2.
- :class:`~starcat.YBSCStarCatalog` for the Yale Bright Star Catalog.
- :class:`~starcat.UCAC4StarCatalog` for UCAC4.

The catalog data itself is not part of this package; see the user guide for where each
catalog looks for its files.
"""

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
