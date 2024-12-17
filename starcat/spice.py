################################################################################
# starcat/spice.py
################################################################################

import os
from pathlib import Path
from typing import Any, Iterator, Optional

import cspyce
from filecache import FCPath

from .starcatalog import Star, StarCatalog


class SpiceStar(Star):
    """A holder for star attributes.

    This class includes attributes unique to stars in SPICE catalogs."""

    def __init__(self) -> None:

        Star.__init__(self)

        self.spectral_class: Optional[str] = None
        """The spectral class"""

        self.temperature: Optional[float] = None
        """The temperature of the star"""


class SpiceStarCatalog(StarCatalog):
    def __init__(self,
                 name: str,
                 dir: Optional[str | Path | FCPath] = None) -> None:

        if dir is None:
            try:
                dir = FCPath(os.environ['SPICE_PATH'])
            except KeyError:
                dir = FCPath(os.environ['OOPS_RESOURCES']) / 'SPICE'
            except KeyError:
                raise RuntimeError(
                    'SPICE_PATH and OOPS_RESOURCES environment variables not set')
        else:
            dir = FCPath(dir)
        self._filename = dir / 'Stars' / f'{name}.bdb'
        local_path = self._filename.retrieve()
        self._catalog = cspyce.stcl01(local_path)[0]
        self.debug_level = 0

        # (ra, dec, ra_uncertainty, dec_uncertainty,
        #  catalog_number, spectral_type, v_magnitude)

    def _find_stars(self,
                    ra_min: float,
                    ra_max: float,
                    dec_min: float,
                    dec_max: float,
                    vmag_min: Optional[float] = None,
                    vmag_max: Optional[float] = None,
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
                star.temperature = self.temperature_from_sclass(star.spectral_class)
                if self.debug_level:
                    print('OK!')

            yield star
