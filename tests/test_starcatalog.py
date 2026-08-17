################################################################################
# tests/test_starcatalog.py
################################################################################

from __future__ import annotations

from typing import Any, Iterator, Optional

import numpy as np
import pytest

from starcat import (SCLASS_TO_B_MINUS_V,
                     SCLASS_TO_SURFACE_TEMP,
                     Star,
                     StarCatalog)
from starcat.starcatalog import (AS_TO_DEG,
                                 AS_TO_RAD,
                                 HALFPI,
                                 MAS_TO_DEG,
                                 MAS_TO_RAD,
                                 TWOPI,
                                 YEAR_TO_SEC)


################################################################################
# Module constants
################################################################################

def test_unit_conversions() -> None:
    assert AS_TO_DEG == pytest.approx(1/3600)
    assert AS_TO_RAD == pytest.approx(np.pi/180/3600)
    assert MAS_TO_DEG == pytest.approx(1/3600/1000)
    assert MAS_TO_RAD == pytest.approx(np.pi/180/3600/1000)
    assert YEAR_TO_SEC == pytest.approx(1/365.25/86400)
    assert TWOPI == pytest.approx(2*np.pi)
    assert HALFPI == pytest.approx(np.pi/2)


def test_sclass_tables() -> None:
    # Both tables must describe the same set of spectral classes
    assert set(SCLASS_TO_B_MINUS_V) == set(SCLASS_TO_SURFACE_TEMP)

    # B-V increases (redder) as the surface temperature decreases
    bmv = list(SCLASS_TO_B_MINUS_V.values())
    temp = list(SCLASS_TO_SURFACE_TEMP.values())
    assert bmv == sorted(bmv)
    assert temp == sorted(temp, reverse=True)


################################################################################
# Star
################################################################################

def test_star_defaults() -> None:
    star = Star()

    assert star.unique_number is None
    assert star.ra is None
    assert star.ra_sigma is None
    assert star.rac_sigma is None
    assert star.dec is None
    assert star.dec_sigma is None
    assert star.vmag is None
    assert star.vmag_sigma is None
    assert star.pm_ra is None
    assert star.pm_ra_sigma is None
    assert star.pm_rac is None
    assert star.pm_rac_sigma is None
    assert star.pm_dec is None
    assert star.pm_dec_sigma is None
    assert star.spectral_class is None
    assert star.temperature is None


def test_star_unknown_attribute() -> None:
    star = Star()

    with pytest.raises(AttributeError) as exc_info:
        _ = star.no_such_attribute

    assert 'Star' in str(exc_info.value)
    assert 'no_such_attribute' in str(exc_info.value)


def test_star_to_from_dict() -> None:
    star = Star()
    star.unique_number = 12345
    star.ra = 1.
    star.dec = -0.5
    star.vmag = 4.25
    star.spectral_class = 'G2'

    d = star.to_dict()
    assert d['unique_number'] == 12345
    assert d['ra'] == 1.
    assert d['dec'] == -0.5
    assert d['vmag'] == 4.25
    assert d['spectral_class'] == 'G2'
    assert d['temperature'] is None
    assert not [key for key in d if key.startswith('__')]
    # Every attribute created by the constructor is present
    assert set(d) == set(vars(Star()))

    star2 = Star()
    star2.from_dict(d)
    assert star2.to_dict() == d


def test_star_from_dict_adds_attributes() -> None:
    star = Star()
    star.from_dict({'unique_number': 7, 'brand_new_attribute': 'hello'})

    assert star.unique_number == 7
    assert star.brand_new_attribute == 'hello'


def test_star_ra_dec_with_pm() -> None:
    star = Star()

    # No RA/DEC at all
    assert star.ra_dec_with_pm(1e6) == (None, None)

    star.ra = 1.
    assert star.ra_dec_with_pm(1e6) == (None, None)

    star.dec = 0.5
    star.ra = None
    assert star.ra_dec_with_pm(1e6) == (None, None)

    # RA/DEC but no proper motion
    star.ra = 1.
    assert star.ra_dec_with_pm(1e6) == (1., 0.5)

    star.pm_ra = 1e-9
    assert star.ra_dec_with_pm(1e6) == (1., 0.5)

    star.pm_ra = None
    star.pm_dec = -2e-9
    assert star.ra_dec_with_pm(1e6) == (1., 0.5)

    # RA/DEC with proper motion
    star.pm_ra = 1e-9
    ra, dec = star.ra_dec_with_pm(1e6)
    assert ra == pytest.approx(1. + 1e-3)
    assert dec == pytest.approx(0.5 - 2e-3)

    # Zero elapsed time is a no-op
    assert star.ra_dec_with_pm(0.) == (1., 0.5)


@pytest.mark.parametrize('sclass,expected',
                         [('A0', 10800),
                          ('a0', 10800),
                          ('  M8  ', 2300),
                          ('O9.5', 31900),
                          ('G2*', 5780),   # SPICE catalogs mark some this way
                          ('Q9', None),
                          ('', None),
                          (None, None)])
def test_temperature_from_sclass(sclass: Optional[str],
                                 expected: Optional[float]) -> None:
    assert Star.temperature_from_sclass(sclass) == expected


@pytest.mark.parametrize('sclass,expected',
                         [('A0', 0.),
                          ('a0', 0.),
                          ('  M8  ', 2.),
                          ('B0.5', -0.28),
                          ('G2*', 0.63),
                          ('Q9', None)])
def test_bmv_from_sclass(sclass: str,
                         expected: Optional[float]) -> None:
    assert Star.bmv_from_sclass(sclass) == expected


def test_sclass_from_bv() -> None:
    # Exact table entries; V is zero so that B-V is exactly B
    assert Star.sclass_from_bv(0., 0.) == 'A0'
    assert Star.sclass_from_bv(0.63, 0.) == 'G2'
    assert Star.sclass_from_bv(-0.32, 0.) == 'O5'   # O5, O6, and O7 tie
    assert Star.sclass_from_bv(2., 0.) == 'M8'

    # Nearest table entry; G2 is +0.63 and G3 is +0.64
    assert Star.sclass_from_bv(0.635, 0.) == 'G2'  # A tie goes to the first
    assert Star.sclass_from_bv(0.6351, 0.) == 'G3'

    # Only the difference matters
    assert Star.sclass_from_bv(10.63, 10.) == 'G2'

    # Outside of the range of the table
    assert Star.sclass_from_bv(-0.33, 0.) is None
    assert Star.sclass_from_bv(2.01, 0.) is None


def test_star_str_full() -> None:
    star = Star()
    star.unique_number = 12345
    star.ra = np.radians(10.5)
    star.ra_sigma = np.radians(1/3600.)
    star.dec = np.radians(20.25)
    star.dec_sigma = np.radians(2/3600.)
    star.vmag = 4.5
    star.vmag_sigma = 0.01
    star.pm_ra = 100. * MAS_TO_RAD * YEAR_TO_SEC
    star.pm_ra_sigma = 5. * MAS_TO_RAD * YEAR_TO_SEC
    star.pm_dec = -200. * MAS_TO_RAD * YEAR_TO_SEC
    star.pm_dec_sigma = 10. * MAS_TO_RAD * YEAR_TO_SEC
    star.temperature = 5780
    star.spectral_class = 'G2'

    s = str(star)

    assert 'UNIQUE ID 12345' in s
    assert 'RA 10.5000000°' in s
    assert '(00h42m0.000s' in s
    assert 'DEC 20.2500000°' in s
    assert '(+020d15m0.000s' in s
    assert '+/- 1.0000s' in s
    assert '+/- 2.0000s' in s
    assert 'VMAG  4.500 +/-  0.010' in s
    assert 'PM RA 100.000 mas/yr +/- 5.000' in s
    assert 'PM DEC -200.000 mas/yr +/- 10.000' in s
    assert 'TEMP  5780' in s
    assert 'SCLASS G2' in s


def test_star_str_minimal() -> None:
    star = Star()
    star.unique_number = 1
    star.spectral_class = 'UNK'

    s = str(star)

    assert 'UNIQUE ID 1' in s
    assert 'RA' not in s
    assert 'DEC' not in s
    assert 'VMAG' not in s
    assert 'PM' not in s
    assert 'TEMP N/A' in s
    assert 'SCLASS UNK' in s


def test_star_str_no_uncertainties() -> None:
    star = Star()
    star.unique_number = 2
    star.ra = np.radians(359.9)
    star.dec = np.radians(-45.5)
    star.vmag = 12.
    star.pm_ra = 0.
    star.pm_dec = 0.
    star.pm_ra_sigma = 0.       # Falsy, so it is not printed
    star.pm_dec_sigma = 0.      # Falsy, so it is not printed
    star.spectral_class = 'K0'

    s = str(star)

    assert 'RA 359.9000000°' in s
    assert 'DEC -45.5000000°' in s
    assert '(-045d30m0.000s)' in s
    assert '+/-' not in s
    assert 'PM RA 0.000 mas/yr' in s
    assert 'PM DEC 0.000 mas/yr' in s
    assert 'TEMP N/A' in s


################################################################################
# StarCatalog
################################################################################

class DummyStarCatalog(StarCatalog):
    """A catalog that records the boxes it is asked to search."""

    def __init__(self, num_stars: int = 1) -> None:
        super().__init__()
        self._num_stars = num_stars
        self.calls: list[tuple[float, float, float, float]] = []
        self.kwargs: list[dict[str, Any]] = []

    def _find_stars(self,
                    ra_min: float,
                    ra_max: float,
                    dec_min: float,
                    dec_max: float,
                    vmag_min: Optional[float] = None,
                    vmag_max: Optional[float] = None,
                    full_result: bool = True,
                    **kwargs: Any) -> Iterator[Star]:

        self.calls.append((ra_min, ra_max, dec_min, dec_max))
        self.kwargs.append({'vmag_min': vmag_min, 'vmag_max': vmag_max,
                            'full_result': full_result, **kwargs})

        for num in range(self._num_stars):
            star = Star()
            star.unique_number = num
            star.ra = ra_min
            star.dec = dec_min
            yield star


def test_star_catalog_defaults() -> None:
    cat = DummyStarCatalog()

    assert cat.debug_level == 0


def test_base_find_stars_not_implemented() -> None:
    cat = StarCatalog()

    with pytest.raises(NotImplementedError):
        list(cat.find_stars())


def test_find_stars_no_split() -> None:
    cat = DummyStarCatalog()

    stars = list(cat.find_stars(ra_min=1., ra_max=2., dec_min=-0.5, dec_max=0.5))

    assert len(stars) == 1
    assert cat.calls == [(1., 2., -0.5, 0.5)]


def test_find_stars_defaults() -> None:
    cat = DummyStarCatalog()

    list(cat.find_stars())

    assert cat.calls == [(0., TWOPI, -HALFPI, HALFPI)]
    assert cat.kwargs == [{'vmag_min': None, 'vmag_max': None,
                           'full_result': True}]


def test_find_stars_clips_out_of_range() -> None:
    cat = DummyStarCatalog()

    list(cat.find_stars(ra_min=-1., ra_max=100., dec_min=-100., dec_max=100.))

    assert cat.calls == [(0., TWOPI, -HALFPI, HALFPI)]


def test_find_stars_split_ra() -> None:
    cat = DummyStarCatalog()

    stars = list(cat.find_stars(ra_min=6., ra_max=1., dec_min=-0.5, dec_max=0.5))

    assert len(stars) == 2
    assert cat.calls == [(0., 1., -0.5, 0.5),
                         (6., TWOPI, -0.5, 0.5)]


def test_find_stars_split_dec() -> None:
    cat = DummyStarCatalog()

    stars = list(cat.find_stars(ra_min=1., ra_max=2., dec_min=0.5, dec_max=-0.5))

    assert len(stars) == 2
    assert cat.calls == [(1., 2., -HALFPI, -0.5),
                         (1., 2., 0.5, HALFPI)]


def test_find_stars_split_ra_and_dec() -> None:
    cat = DummyStarCatalog()

    stars = list(cat.find_stars(ra_min=6., ra_max=1., dec_min=0.5, dec_max=-0.5))

    assert len(stars) == 4
    assert cat.calls == [(0., 1., -HALFPI, -0.5),
                         (6., TWOPI, -HALFPI, -0.5),
                         (0., 1., 0.5, HALFPI),
                         (6., TWOPI, 0.5, HALFPI)]


@pytest.mark.parametrize('ra_min,ra_max,dec_min,dec_max,num_calls',
                         [(1., 2., -0.5, 0.5, 1),
                          (6., 1., -0.5, 0.5, 2),
                          (1., 2., 0.5, -0.5, 2),
                          (6., 1., 0.5, -0.5, 4)])
def test_find_stars_passes_arguments(ra_min: float,
                                     ra_max: float,
                                     dec_min: float,
                                     dec_max: float,
                                     num_calls: int) -> None:
    cat = DummyStarCatalog()

    list(cat.find_stars(ra_min=ra_min, ra_max=ra_max,
                        dec_min=dec_min, dec_max=dec_max,
                        vmag_min=1., vmag_max=10., full_result=False,
                        extra_option='surprise'))

    assert len(cat.kwargs) == num_calls
    for kwargs in cat.kwargs:
        assert kwargs == {'vmag_min': 1., 'vmag_max': 10.,
                          'full_result': False, 'extra_option': 'surprise'}


def test_count_stars() -> None:
    cat = DummyStarCatalog(num_stars=3)

    assert cat.count_stars() == 3
    assert cat.kwargs == [{'vmag_min': None, 'vmag_max': None,
                           'full_result': False}]

    # Four boxes searched, three stars in each
    cat = DummyStarCatalog(num_stars=3)
    assert cat.count_stars(ra_min=6., ra_max=1., dec_min=0.5, dec_max=-0.5) == 12


def test_count_stars_empty() -> None:
    cat = DummyStarCatalog(num_stars=0)

    assert cat.count_stars() == 0
