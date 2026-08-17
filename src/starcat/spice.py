################################################################################
# starcat/spice.py
################################################################################

"""Support for NAIF SPICE type 1 star catalogs.

A SPICE star catalog is a single ``.bdb`` kernel read through the CSPICE ``STC01``
routines; catalogs distributed by NAIF include ``hipparcos``, ``ppm``, and ``tycho2``.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import cspyce
from filecache import FCPath

from .starcatalog import Star, StarCatalog


class SpiceStar(Star):
    """A holder for star attributes.

    This class includes attributes unique to stars in SPICE catalogs.

    A SpiceStar only supports these attributes: `unique_number`, `ra`, `ra_sigma`, `dec`,
    `dec_sigma`, `vmag`, `spectral_class`, `temperature`
    """

    def __init__(self) -> None:
        """Constructor for SpiceStar."""

        # Initialize the standard fields
        super().__init__()


class SpiceStarCatalog(StarCatalog):
    """A NAIF SPICE type 1 star catalog.

    The catalog is a single ``.bdb`` kernel, which is loaded when this object is created
    and stays loaded for the life of the process. `find_stars` and `count_stars` accept no
    options beyond the ones common to all catalogs, and the resulting :class:`SpiceStar`
    objects carry a catalog number, a position and its uncertainty, a magnitude, a spectral
    class, and a temperature derived from that class. SPICE catalogs record no proper
    motion.
    """

    def __init__(self,
                 name: str,
                 dir: str | Path | FCPath | None = None) -> None:
        """Create a SpiceStarCatalog.

        Parameters:
            name: The name of the SPICE catalog without the extension, such as
                ``hipparcos``, ``ppm``, or ``tycho2``.
            dir: The path to the star catalog directory (may be a URL). Within this
                directory should be the kernel for the requested name (``name.bdb``). If
                None, the ``Stars`` subdirectory of ``SPICE_PATH`` is used, or failing
                that the ``SPICE/Stars`` subdirectory of ``OOPS_RESOURCES``.

        Raises:
            RuntimeError: If ``dir`` is None and neither environment variable is set.
        """

        super().__init__()
        if dir is None:
            if 'SPICE_PATH' in os.environ:
                dir = FCPath(os.environ['SPICE_PATH']) / 'Stars'
            elif 'OOPS_RESOURCES' in os.environ:
                dir = FCPath(os.environ['OOPS_RESOURCES']) / 'SPICE' / 'Stars'
            else:
                raise RuntimeError(
                    'SPICE_PATH and OOPS_RESOURCES environment variables not set')
        else:
            dir = FCPath(dir)
        self._filename = dir / f'{name}.bdb'
        local_path = self._filename.retrieve()
        self._catalog = cspyce.stcl01(local_path)[0]

    def _find_stars(self,
                    ra_min: float,
                    ra_max: float,
                    dec_min: float,
                    dec_max: float,
                    vmag_min: float | None = None,
                    vmag_max: float | None = None,
                    full_result: bool = True,
                    **kwargs: Any) -> Iterator[SpiceStar]:

        nstars = cspyce.stcf01(self._catalog, ra_min, ra_max, dec_min, dec_max)

        for i in range(nstars):
            star = SpiceStar()
            result = tuple(cspyce.stcg01(i))
            (star.ra, star.dec, star.ra_sigma, star.dec_sigma,
             star.unique_number, star.spectral_class, star.vmag) = result
            if star.vmag is not None:
                if vmag_min is not None and star.vmag < vmag_min:
                    if self.debug_level:
                        print('SKIPPED VMAG', star.vmag)
                    continue
                if vmag_max is not None and star.vmag > vmag_max:
                    if self.debug_level:
                        print('SKIPPED VMAG', star.vmag)
                    continue

            if full_result:
                star.temperature = Star.temperature_from_sclass(star.spectral_class)

            if self.debug_level:
                print(star)
                print('-' * 80)

            yield star
