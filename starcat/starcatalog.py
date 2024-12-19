################################################################################
# starcat/starcatalog.py
################################################################################

from __future__ import annotations

import inspect
import numpy as np
from typing import Any, Iterator, Optional


AS_TO_DEG = 1 / 3600.
AS_TO_RAD = np.radians(AS_TO_DEG)
MAS_TO_DEG = AS_TO_DEG / 1000.
MAS_TO_RAD = np.radians(MAS_TO_DEG)
YEAR_TO_SEC = 1 / 365.25 / 86400.

TWOPI = 2 * np.pi
HALFPI = np.pi / 2

#===============================================================================
#
# Jacobson B-V photometry vs. stellar spectral classification
#
# Data from
#     Zombeck, M. V. Handbook of Space Astronomy and Astrophysics,
#     Cambridge, UK: Cambridge University Press, 2nd ed., pp. 68-70
#     http://ads.harvard.edu/books/hsaa/
# as transcribed:
#     http://www.vendian.org/mncharity/dir3/starcolor/details.html
#
# The tables below only include main sequence stars, but the temperature
# difference between main sequence stars and giant stars is minimal for our
# purposes. Missing values have been linearly interpolated.
#
#===============================================================================

SCLASS_TO_B_MINUS_V = {
    'O5': -0.32,
    'O6': -0.32,
    'O7': -0.32,
    'O8': -0.31,
    'O9': -0.31,
    'O9.5': -0.30,
    'B0': -0.30,
    'B0.5': -0.28,
    'B1': -0.26,
    'B2': -0.24,
    'B3': -0.20,
    'B4': -0.18,  # Interpolated
    'B5': -0.16,
    'B6': -0.14,
    'B7': -0.12,
    'B8': -0.09,
    'B9': -0.06,
    'A0': +0.00,
    'A1': +0.02,  # Interpolated
    'A2': +0.04,  # Interpolated
    'A3': +0.06,
    'A4': +0.10,  # Interpolated
    'A5': +0.14,
    'A6': +0.16,  # Interpolated
    'A7': +0.19,
    'A8': +0.23,  # Interpolated
    'A9': +0.27,  # Interpolated
    'F0': +0.31,
    'F1': +0.33,  # Interpolated
    'F2': +0.36,
    'F3': +0.38,  # Interpolated
    'F4': +0.41,  # Interpolated
    'F5': +0.43,
    'F6': +0.47,  # Interpolated
    'F7': +0.51,  # Interpolated
    'F8': +0.54,
    'F9': +0.56,  # Interpolated
    'G0': +0.59,
    'G1': +0.61,  # Interpolated
    'G2': +0.63,
    'G3': +0.64,  # Interpolated
    'G4': +0.65,  # Interpolated
    'G5': +0.66,
    'G6': +0.69,  # Interpolated
    'G7': +0.72,  # Interpolated
    'G8': +0.74,
    'G9': +0.78,  # Interpolated
    'K0': +0.82,
    'K1': +0.87,  # Interpolated
    'K2': +0.92,
    'K3': +0.99,  # Interpolated
    'K4': +1.07,  # Interpolated
    'K5': +1.15,
    'K6': +1.22,  # Interpolated
    'K7': +1.30,
    'K8': +1.33,  # Interpolated
    'K9': +1.37,  # Interpolated
    'M0': +1.41,
    'M1': +1.48,
    'M2': +1.52,
    'M3': +1.55,
    'M4': +1.56,
    'M5': +1.61,
    'M6': +1.72,
    'M7': +1.84,
    'M8': +2.00
}

SCLASS_TO_SURFACE_TEMP = {
    'O5': 38000,
    'O6': 38000,
    'O7': 38000,
    'O8': 35000,
    'O9': 35000,
    'O9.5': 31900,
    'B0': 30000,
    'B0.5': 27000,
    'B1': 24200,
    'B2': 22100,
    'B3': 18800,
    'B4': 17600,  # Interpolated
    'B5': 16400,
    'B6': 15400,
    'B7': 14500,
    'B8': 13400,
    'B9': 12400,
    'A0': 10800,
    'A1': 10443,  # Interpolated
    'A2': 10086,  # Interpolated
    'A3': 9730,
    'A4': 9175,  # Interpolated
    'A5': 8620,
    'A6': 8405,  # Interpolated
    'A7': 8190,
    'A8': 7873,  # Interpolated
    'A9': 7557,  # Interpolated
    'F0': 7240,
    'F1': 7085,  # Interpolated
    'F2': 6930,
    'F3': 6800,  # Interpolated
    'F4': 6670,  # Interpolated
    'F5': 6540,
    'F6': 6427,  # Interpolated
    'F7': 6313,  # Interpolated
    'F8': 6200,
    'F9': 6060,  # Interpolated
    'G0': 5920,
    'G1': 5850,  # Interpolated
    'G2': 5780,
    'G3': 5723,  # Interpolated
    'G4': 5667,  # Interpolated
    'G5': 5610,
    'G6': 5570,  # Interpolated
    'G7': 5530,  # Interpolated
    'G8': 5490,
    'G9': 5365,  # Interpolated
    'K0': 5240,
    'K1': 5010,  # Interpolated
    'K2': 4780,
    'K3': 4706,  # Interpolated
    'K4': 4632,  # Interpolated
    'K5': 4558,  # Interpolated
    'K6': 4484,  # Interpolated
    'K7': 4410,
    'K8': 4247,  # Interpolated
    'K9': 4083,  # Interpolated
    'M0': 3800,  # M class from https://arxiv.org/abs/0903.3371
    'M1': 3600,
    'M2': 3400,
    'M3': 3250,
    'M4': 3100,
    'M5': 2800,
    'M6': 2600,
    'M7': 2500,
    'M8': 2300,
}


#===============================================================================
#
# STAR Superclass
#
#===============================================================================

class Star:
    """A holder for star attributes.

    This is the base class that defines attributes common to all or most
    star catalogs. Note that fields vary among star catalogs and individual
    stars and there is no guarantee a particular field will be filled in.
    """

    def __init__(self) -> None:
        """Constructor for Star; additional attributes are available from subclasses."""

        self.unique_number: Optional[int] = None
        """Unique catalog number (may not be unique across catalogs)"""

        self.ra: Optional[float] = None
        """Right ascension at J2000 epoch (radians)"""

        self.ra_sigma: Optional[float] = None
        """Right ascension uncertainty (radians)"""

        self.rac_sigma: Optional[float] = None
        """Right ascension * cos(DEC) uncertainty (radians)"""

        self.dec: Optional[float] = None
        """Declination at J2000 epoch (radians)"""

        self.dec_sigma: Optional[float] = None
        """Declination uncertainty (radians)"""

        self.vmag: Optional[float] = None
        """Visual magnitude"""

        self.vmag_sigma: Optional[float] = None
        """Visual magnitude uncertainty"""

        self.pm_ra: Optional[float] = None
        """Proper motion in RA (radians/sec)"""

        self.pm_ra_sigma: Optional[float] = None
        """Proper motion in RA uncertainty (radians/sec)"""

        self.pm_rac: Optional[float] = None
        """Proper motion in RA * cos(DEC) (radians/sec)"""

        self.pm_rac_sigma: Optional[float] = None
        """Proper motion in RA * cos(DEC) uncertainty (radians/sec)"""

        self.pm_dec: Optional[float] = None
        """Proper motion in DEC (radians/sec)"""

        self.pm_dec_sigma: Optional[float] = None
        """Proper motion in DEC error (radians/sec)"""

        self.spectral_class: Optional[str] = None
        """Spectral class"""

        self.temperature: Optional[float] = None
        """Star temperature (usually derived from spectral class)"""

    def __str__(self) -> str:

        ret = f'UNIQUE ID {self.unique_number:d}'

        if self.ra is not None:
            ret += f' | RA {np.degrees(self.ra):.7f}°'
            if self.ra_sigma is not None:
                ret += f' [+/- {np.degrees(self.ra_sigma):.7f}°]'

            ra_deg = np.degrees(self.ra)/15  # In hours
            hh = int(ra_deg)
            mm = int((ra_deg-hh)*60)
            ss = (ra_deg-hh-mm/60.)*3600
            ret += f' ({hh:02d}h{mm:02d}m{ss:05.3f}s'
            if self.ra_sigma is not None:
                ret += ' +/- %.4fs' % (np.degrees(self.ra_sigma)*3600)
            ret += ')'

        if self.dec is not None:
            ret += f' | DEC {np.degrees(self.dec):.7f}°'
            if self.dec_sigma is not None:
                ret += f' [+/- {np.degrees(self.dec_sigma):.7f}°]'

            dec_deg = np.degrees(self.dec)
            neg = '+'
            if dec_deg < 0.:
                neg = '-'
                dec_deg = -dec_deg
            dd = int(dec_deg)
            mm = int((dec_deg-dd)*60)
            ss = (dec_deg-dd-mm/60.)*3600
            ret += f' ({neg}{dd:03d}d{mm:02d}m{ss:05.3f}s'

            if self.dec_sigma is not None:
                ret += ' +/- %.4fs' % (np.degrees(self.dec_sigma)*3600)
            ret += ')'

        ret += '\n'

        if self.vmag is not None:
            ret += f'VMAG {self.vmag:6.3f} '
            if self.vmag_sigma is not None:
                ret += f'+/- {self.vmag_sigma:6.3f} '

        if self.pm_ra is not None:
            ret += ' | PM RA %.3f mas/yr ' % (self.pm_ra / MAS_TO_RAD / YEAR_TO_SEC)
            if self.pm_ra_sigma:
                ret += '+/- %.3f ' % (self.pm_ra_sigma / MAS_TO_RAD / YEAR_TO_SEC)

        if self.pm_dec is not None:
            ret += ' | PM DEC %.3f mas/yr ' % (self.pm_dec / MAS_TO_RAD / YEAR_TO_SEC)
            if self.pm_dec_sigma:
                ret += '+/- %.3f ' % (self.pm_dec_sigma / MAS_TO_RAD / YEAR_TO_SEC)

        ret += '\n'

        if self.temperature is None:
            ret += 'TEMP N/A'
        else:
            ret += f'TEMP {self.temperature:5d}'
        ret += f' | SCLASS {self.spectral_class:s}'

        return ret

    def to_dict(self) -> dict[str, Any]:
        """Return a dictionary containing all star attributes."""

        attribs = inspect.getmembers(self, lambda a: not inspect.isroutine(a))
        attribs = [a for a in attribs
                   if not (a[0].startswith('__') and a[0].endswith('__'))]
        return dict(attribs)

    def from_dict(self, d: dict[str, Any]) -> None:
        """Set the attributes for this star based on the dictionary."""

        for key in list(d.keys()):
            setattr(self, key, d[key])

    def ra_dec_with_pm(self, tdb: float) -> tuple[float, float] | tuple[None, None]:
        """Return the star's RA and DEC adjusted for proper motion.

        If no proper motion is available, the original RA and DEC are returned.

        Parameters:
            tdb: The time since the J2000 epoch in seconds.

        Returns:
            A tuple containing the RA and DEC adjusted for proper motion, if possible.
        """

        if self.ra is None or self.dec is None:
            return (None, None)

        if self.pm_ra is None or self.pm_dec is None:
            return (self.ra, self.dec)

        return (self.ra + tdb*self.pm_ra, self.dec + tdb*self.pm_dec)

    @staticmethod
    def temperature_from_sclass(sclass: Optional[str]) -> float | None:
        """Return a star's temperature (K) given its spectral class.

        Parameters:
            sclass: The spectral class.

        Returns:
            The temperature associated with the spectral class.
        """

        if sclass is None:
            return None
        if sclass.endswith('*'):  # This happens on some SPICE catalog stars
            sclass = sclass[:-1]
        sclass = sclass.strip().upper()
        try:
            return SCLASS_TO_SURFACE_TEMP[sclass]
        except KeyError:
            return None

    @staticmethod
    def bmv_from_sclass(sclass: str) -> float | None:
        """Return a star's B-V color given its spectral class.

        Parameters:
            sclass: The spectral class.

        Returns:
            The B-V color associated with the spectral class.
        """

        if sclass[-1] == '*':  # This happens on some SPICE catalog stars
            sclass = sclass[:-1]
        sclass = sclass.strip().upper()
        try:
            return SCLASS_TO_B_MINUS_V[sclass]
        except KeyError:
            return None

    @staticmethod
    def sclass_from_bv(b: float,
                       v: float) -> str | None:
        """Return a star's spectral class given photometric B and V.

        Parameters:
            b: The photometric B value.
            v: The photometric V value.

        Returns:
            The spectral class, if available.
        """

        bmv = b - v

        best_sclass = None
        best_resid = 1e38

        min_bmv = 1e38
        max_bmv = -1e38
        for sclass, sbmv in SCLASS_TO_B_MINUS_V.items():
            min_bmv = min(min_bmv, sbmv)
            max_bmv = max(max_bmv, sbmv)
            resid = abs(sbmv-bmv)
            if resid < best_resid:
                best_resid = resid
                best_sclass = sclass

        if min_bmv <= bmv <= max_bmv:
            return best_sclass

        return None


class StarCatalog:
    def __init__(self) -> None:
        self.debug_level = 0

    def count_stars(self, **kwargs: Any) -> int:
        """Count the stars that match the given search criteria."""

        count = 0
        for _ in self.find_stars(full_result=False, **kwargs):
            count += 1
        return count

    def find_stars(self,
                   ra_min: float = 0,
                   ra_max: float = TWOPI,
                   dec_min: float = -HALFPI,
                   dec_max: float = HALFPI,
                   vmag_min: Optional[float] = None,
                   vmag_max: Optional[float] = None,
                   full_result: bool = True,
                   **kwargs: Any) -> Iterator[Star]:
        """Yield the stars that match the given search criteria.

        Parameters:
            ra_min: The minimum RA.
            ra_max: The maximum RA.
            dec_min: The minimum DEC.
            dec_max: The maximum DEC.
            vmag_min: The minimum visual magnitude.
            vmag_max: The maximum visual magnitude.
            full_result: If True, fill in all available fields of the resulting
                :class:`Star`. If False, some fields will not be filled in to save
                time. This is most useful when counting stars.

        Yields:
            The :class:`Star` objects that meet the given constraints.
        """

        ra_min = np.clip(ra_min, 0., TWOPI)
        ra_max = np.clip(ra_max, 0., TWOPI)
        dec_min = np.clip(dec_min, -HALFPI, HALFPI)
        dec_max = np.clip(dec_max, -HALFPI, HALFPI)

        if ra_min > ra_max:
            if dec_min > dec_max:
                # Split into four searches
                for star in self._find_stars(0., ra_max, -HALFPI, dec_max,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
                for star in self._find_stars(ra_min, TWOPI, -HALFPI, dec_max,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
                for star in self._find_stars(0., ra_max, dec_min, HALFPI,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
                for star in self._find_stars(ra_min, TWOPI, dec_min, HALFPI,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
            else:
                # Split into two searches - RA
                for star in self._find_stars(0., ra_max, dec_min, dec_max,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
                for star in self._find_stars(ra_min, TWOPI, dec_min, dec_max,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
        else:
            if dec_min > dec_max:
                # Split into two searches - DEC
                for star in self._find_stars(ra_min, ra_max, -HALFPI, dec_max,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
                for star in self._find_stars(ra_min, ra_max, dec_min, HALFPI,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star
            else:
                # No need to split at all
                for star in self._find_stars(ra_min, ra_max, dec_min, dec_max,
                                             vmag_min=vmag_min, vmag_max=vmag_max,
                                             full_result=full_result,
                                             **kwargs):
                    yield star

    def _find_stars(self,
                    ra_min: float,
                    ra_max: float,
                    dec_min: float,
                    dec_max: float,
                    vmag_min: Optional[float] = None,
                    vmag_max: Optional[float] = None,
                    full_result: bool = True,
                    **kwargs: Any) -> Iterator[Star]:
        raise NotImplementedError
