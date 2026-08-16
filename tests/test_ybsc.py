################################################################################
# tests/test_ybsc.py
################################################################################

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from catalog_data import write_ybsc_catalog, ybsc_dec, ybsc_ra, ybsc_record
from filecache import FCPath

from starcat import Star, YBSCStar, YBSCStarCatalog
from starcat.starcatalog import AS_TO_RAD, YEAR_TO_SEC

# A Vega-like star with every optional field filled in
VEGA = ybsc_record(hr='7001', name='3Alp Lyr', dm='BD+38 3238', hd='172167',
                   sao='67174', fk5='699', ir_flag='I', ir_ref=' ',
                   multiple=' ', ads='11510', ads_comp='AB', var_id='Alp Lyr',
                   **ybsc_ra(18, 36, 56.3), **ybsc_dec('+', 38, 47, 1),
                   glon='67.44', glat='19.24',
                   vmag='0.03', vmag_code=' ', vmag_uncertainty=' ',
                   b_v='0.00', u_b='-0.01', r_i='-0.03', r_i_code='C',
                   sptype='A0Va', sptype_code='v',
                   pm_ra='0.202', pm_dec='0.286',
                   parallax_type='D', parallax='0.123',
                   rad_vel='-014', rad_vel_comments='V',
                   rot_vel_limit='<', rot_vel='15', rot_vel_uncertainty='v',
                   dmag='10.4', sep='62.8', mult_id='AB', mult_cnt='5',
                   note_flag='*')

VEGA_RA = np.radians((18 + 36/60. + 56.3/3600.) * 15.)
VEGA_DEC = np.radians(38 + 47/60. + 1/3600.)


def make_catalog(root: Path,
                 records: Sequence[str]) -> YBSCStarCatalog:
    """Write a YBSC catalog and return a catalog object for it."""

    write_ybsc_catalog(root, records)

    return YBSCStarCatalog(root)


def around(cat: YBSCStarCatalog,
           ra: float,
           dec: float,
           **kwargs: Any) -> list[Star]:
    """Return the stars within 0.1 degrees of the given RA and DEC."""

    delta = np.radians(0.1)

    return list(cat.find_stars(ra_min=ra-delta, ra_max=ra+delta,
                               dec_min=dec-delta, dec_max=dec+delta,
                               **kwargs))


################################################################################
# YBSCStar
################################################################################

def test_ybsc_star_defaults() -> None:
    star = YBSCStar()

    # Inherited from Star
    assert star.unique_number is None
    assert star.ra is None
    assert star.dec is None

    # YBSC-specific
    assert star.name is None
    assert star.durchmusterung_id is None
    assert star.draper_number is None
    assert star.sao_number is None
    assert star.fk5_number is None
    assert star.ir_source is None
    assert star.ir_source_ref is None
    assert star.multiple_star_code is None
    assert star.aitken_designation is None
    assert star.ads_components is None
    assert star.variable_star_id is None
    assert star.galactic_longitude is None
    assert star.galactic_latitude is None
    assert star.vmag_code is None
    assert star.vmag_uncertainty_flag is None
    assert star.b_v is None
    assert star.b_v_uncertainty_flag is None
    assert star.u_b is None
    assert star.u_b_uncertainty_flag is None
    assert star.r_i is None
    assert star.r_i_code is None
    assert star.spectral_class_code is None
    assert star.parallax_type is None
    assert star.parallax is None
    assert star.radial_velocity is None
    assert star.radial_velocity_comments is None
    assert star.rotational_velocity_limit is None
    assert star.rotational_velocity is None
    assert star.rotational_velocity_uncertainty_flag is None
    assert star.double_mag_diff is None
    assert star.double_mag_sep is None
    assert star.double_mag_components is None
    assert star.multiple_num_components is None


################################################################################
# Reading records
################################################################################

def test_parse_star(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [VEGA])

    assert len(cat._stars) == 1
    star = cat._stars[0]

    assert isinstance(star, YBSCStar)

    # Catalog numbers
    assert star.unique_number == 7001
    assert star.name == '3Alp Lyr'
    assert star.durchmusterung_id == 'BD+38 3238'
    assert star.draper_number == 172167
    assert star.sao_number == 67174
    assert star.fk5_number == 699

    # Source flags
    assert star.ir_source is True
    assert star.ir_source_ref == YBSCStar.YBSC_IR_NASA
    assert star.multiple_star_code == ' '
    assert star.aitken_designation == '11510'
    assert star.ads_components == 'AB'
    assert star.variable_star_id == 'Alp Lyr'

    # Position
    assert star.ra == pytest.approx(VEGA_RA)
    assert star.dec == pytest.approx(VEGA_DEC)
    assert star.galactic_longitude == pytest.approx(np.radians(67.44))
    assert star.galactic_latitude == pytest.approx(np.radians(19.24))

    # Magnitudes
    assert star.vmag == pytest.approx(0.03)
    assert star.vmag_code == ' '
    assert star.vmag_uncertainty_flag == ' '
    assert star.b_v == pytest.approx(0.)
    assert star.u_b == pytest.approx(-0.01)
    assert star.r_i == pytest.approx(-0.03)
    assert star.r_i_code == 'C'

    # Spectral class
    assert star.spectral_class == 'A0Va'
    assert star.spectral_class_code == 'v'
    assert star.temperature == 10800

    # Motion and parallax
    assert star.pm_rac == pytest.approx(0.202 * AS_TO_RAD * YEAR_TO_SEC)
    assert star.pm_ra == pytest.approx(0.202 * AS_TO_RAD * YEAR_TO_SEC /
                                       np.cos(VEGA_DEC))
    assert star.pm_dec == pytest.approx(0.286 * AS_TO_RAD * YEAR_TO_SEC)
    assert star.parallax_type == 'D'
    assert star.parallax == pytest.approx(0.123)
    assert star.radial_velocity == pytest.approx(-14.)
    assert star.radial_velocity_comments == 'V'
    assert star.rotational_velocity_limit == '<'
    assert star.rotational_velocity == pytest.approx(15.)
    assert star.rotational_velocity_uncertainty_flag == 'v'

    # Doubles and multiples
    assert star.double_mag_diff == pytest.approx(10.4)
    assert star.double_mag_sep == pytest.approx(62.8 * AS_TO_RAD)
    assert star.double_mag_components == 'AB'
    assert star.multiple_num_components == 5


def test_parse_minimal_star(tmp_path: Path) -> None:
    # Only the fields that every record must have
    minimal = ybsc_record(hr='1', **ybsc_ra(0, 0, 1.1),
                          **ybsc_dec('+', 44, 40, 22),
                          glon='114.44', glat='-16.88', vmag='6.70',
                          sptype='A1Vn', pm_ra='-0.012', pm_dec='-0.018')
    cat = make_catalog(tmp_path, [minimal])

    star = cat._stars[0]

    assert star.unique_number == 1
    assert star.name == ''
    assert star.durchmusterung_id == ''
    assert star.draper_number is None
    assert star.sao_number is None
    assert star.fk5_number is None
    assert star.ir_source is False
    assert star.ir_source_ref == YBSCStar.YBSC_IR_NASA
    assert star.aitken_designation is None
    assert star.ads_components is None
    assert star.variable_star_id is None
    assert star.b_v is None
    assert star.u_b is None
    assert star.r_i is None
    assert star.parallax is None
    assert star.radial_velocity is None
    assert star.radial_velocity_comments == ''
    assert star.rotational_velocity is None
    assert star.rotational_velocity_limit == ''
    assert star.double_mag_diff is None
    assert star.double_mag_sep is None
    assert star.double_mag_components == ''
    assert star.multiple_num_components is None

    assert star.vmag == pytest.approx(6.70)
    assert star.spectral_class == 'A1Vn'
    assert star.temperature == 10443


def test_record_shorter_than_full_width(tmp_path: Path) -> None:
    # Records are padded on the right, as the real catalog is ragged
    short = ybsc_record(hr='1', **ybsc_ra(0, 0, 1.1),
                        **ybsc_dec('+', 44, 40, 22),
                        glon='114.44', glat='-16.88', vmag='6.70',
                        sptype='A1Vn', pm_ra='-0.012', pm_dec='-0.018')
    short = short.rstrip()
    assert len(short) < 197

    cat = make_catalog(tmp_path, [short])

    assert len(cat._stars) == 1
    assert cat._stars[0].multiple_num_components is None


def test_records_without_vmag_are_skipped(tmp_path: Path) -> None:
    # Stars removed from the catalog have blank position and magnitude fields
    removed = ybsc_record(hr='92')
    cat = make_catalog(tmp_path, [VEGA, removed])

    assert len(cat._stars) == 1
    assert cat._stars[0].unique_number == 7001


@pytest.mark.parametrize(('flag', 'expected'),
                         [(' ', YBSCStar.YBSC_IR_NASA),
                          ("'", YBSCStar.YBSC_IR_ENGLES),
                          (':', YBSCStar.YBSC_IR_UNCERTAIN),
                          ('X', None)])
def test_ir_source_reference(tmp_path: Path,
                             flag: str,
                             expected: int | None) -> None:
    cat = make_catalog(tmp_path, [ybsc_record(hr='1', ir_flag='I', ir_ref=flag,
                                              **ybsc_ra(1, 0, 0.),
                                              **ybsc_dec('+', 10, 0, 0),
                                              glon='0.00', glat='0.00',
                                              vmag='5.00', sptype='G0',
                                              pm_ra='0.000', pm_dec='0.000')])

    star = cat._stars[0]

    assert star.ir_source is True
    assert star.ir_source_ref == expected


@pytest.mark.parametrize(('sptype', 'temperature'),
                         [('A0Va', 10800),
                          ('gG9', 5365),      # A giant; the g is stripped
                          ('K0IIIbCN-0.5', 5240),
                          ('QQQ', None)])
def test_temperature_from_spectral_type(tmp_path: Path,
                                        sptype: str,
                                        temperature: int | None) -> None:
    cat = make_catalog(tmp_path, [ybsc_record(hr='1', **ybsc_ra(1, 0, 0.),
                                              **ybsc_dec('+', 10, 0, 0),
                                              glon='0.00', glat='0.00',
                                              vmag='5.00', sptype=sptype,
                                              pm_ra='0.000', pm_dec='0.000')])

    star = cat._stars[0]

    assert star.spectral_class == sptype
    assert star.temperature == temperature


@pytest.mark.parametrize(('sign', 'degrees', 'minutes', 'seconds', 'expected'),
                         [('+', 38, 47, 1, 38 + 47/60. + 1/3600.),
                          ('-', 38, 47, 1, -(38 + 47/60. + 1/3600.)),
                          ('+', 0, 30, 11, 30/60. + 11/3600.),
                          ('-', 0, 30, 11, -(30/60. + 11/3600.)),
                          ('+', 89, 59, 59, 89 + 59/60. + 59/3600.),
                          ('-', 89, 59, 59, -(89 + 59/60. + 59/3600.))])
def test_declination_sign(tmp_path: Path,
                          sign: str,
                          degrees: int,
                          minutes: int,
                          seconds: int,
                          expected: float) -> None:
    # A star less than one degree south of the equator has a degrees field of
    # zero, so only the sign character distinguishes it from a northern star
    cat = make_catalog(tmp_path,
                       [ybsc_record(hr='1', **ybsc_ra(1, 0, 0.),
                                    **ybsc_dec(sign, degrees, minutes, seconds),
                                    glon='0.00', glat='0.00', vmag='5.00',
                                    sptype='G0', pm_ra='0.000',
                                    pm_dec='0.000')])

    assert cat._stars[0].dec == pytest.approx(np.radians(expected))


def test_right_ascension(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [ybsc_record(hr='1', **ybsc_ra(23, 59, 59.9),
                                              **ybsc_dec('+', 10, 0, 0),
                                              glon='0.00', glat='0.00',
                                              vmag='5.00', sptype='G0',
                                              pm_ra='0.000', pm_dec='0.000')])

    assert cat._stars[0].ra == pytest.approx(
        np.radians((23 + 59/60. + 59.9/3600.) * 15.))


################################################################################
# Searching
################################################################################

def test_find_stars_by_position(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [VEGA])

    assert len(around(cat, VEGA_RA, VEGA_DEC)) == 1
    assert around(cat, VEGA_RA, VEGA_DEC)[0].vmag == pytest.approx(0.03)

    # Just outside the box in each direction
    offset = np.radians(0.3)
    assert around(cat, VEGA_RA+offset, VEGA_DEC) == []
    assert around(cat, VEGA_RA-offset, VEGA_DEC) == []
    assert around(cat, VEGA_RA, VEGA_DEC+offset) == []
    assert around(cat, VEGA_RA, VEGA_DEC-offset) == []


def test_find_stars_by_vmag(tmp_path: Path) -> None:
    stars = [ybsc_record(hr=str(num+1), **ybsc_ra(1, 0, 0.),
                         **ybsc_dec('+', 10, 0, 0), glon='0.00', glat='0.00',
                         vmag=f'{vmag:5.2f}', sptype='G0',
                         pm_ra='0.000', pm_dec='0.000')
             for num, vmag in enumerate([1., 5., 9.])]
    cat = make_catalog(tmp_path, stars)

    assert cat.count_stars() == 3
    assert cat.count_stars(vmag_min=3.) == 2
    assert cat.count_stars(vmag_max=7.) == 2
    assert cat.count_stars(vmag_min=3., vmag_max=7.) == 1


def test_find_stars_by_vmag_of_zero(tmp_path: Path) -> None:
    # A limit of zero is a real limit, not the absence of one; the bright stars of
    # this catalog straddle magnitude zero
    stars = [ybsc_record(hr=str(num+1), **ybsc_ra(1, 0, 0.),
                         **ybsc_dec('+', 10, 0, 0), glon='0.00', glat='0.00',
                         vmag=f'{vmag:5.2f}', sptype='G0',
                         pm_ra='0.000', pm_dec='0.000')
             for num, vmag in enumerate([-1., 1.])]
    cat = make_catalog(tmp_path, stars)

    assert cat.count_stars() == 2
    assert cat.count_stars(vmag_min=0.) == 1
    assert cat.count_stars(vmag_max=0.) == 1


def test_allow_double(tmp_path: Path) -> None:
    single = ybsc_record(hr='1', **ybsc_ra(1, 0, 0.),
                         **ybsc_dec('+', 10, 0, 0), glon='0.00', glat='0.00',
                         vmag='5.00', sptype='G0',
                         pm_ra='0.000', pm_dec='0.000')
    double = ybsc_record(hr='2', multiple='W', **ybsc_ra(1, 0, 0.),
                         **ybsc_dec('+', 10, 0, 0), glon='0.00', glat='0.00',
                         vmag='5.00', sptype='G0',
                         pm_ra='0.000', pm_dec='0.000')
    cat = make_catalog(tmp_path, [single, double])

    assert len(cat._stars) == 2
    assert cat.count_stars() == 1
    assert cat.count_stars(allow_double=True) == 2
    assert [s.unique_number for s in cat.find_stars(allow_double=True)] == [1, 2]


def test_stars_without_position_or_magnitude(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [VEGA])

    no_ra = YBSCStar()
    no_ra.dec = 0.
    no_ra.multiple_star_code = ' '
    no_dec = YBSCStar()
    no_dec.ra = 0.
    no_dec.multiple_star_code = ' '
    no_vmag = YBSCStar()
    no_vmag.ra = VEGA_RA
    no_vmag.dec = VEGA_DEC
    no_vmag.multiple_star_code = ' '
    cat._stars.extend([no_ra, no_dec, no_vmag])

    # The stars with no position are ignored; the one with no magnitude is
    # returned no matter what magnitude limits are given
    assert cat.count_stars() == 2
    assert cat.count_stars(vmag_min=3.) == 1
    assert cat.count_stars(vmag_max=3.) == 2


def test_debug_output(tmp_path: Path,
                      capsys: pytest.CaptureFixture[str]) -> None:
    cat = make_catalog(tmp_path, [VEGA])
    cat.debug_level = 1

    assert cat.count_stars() == 1

    out = capsys.readouterr().out
    assert 'UNIQUE ID 7001' in out
    assert '-' * 80 in out


################################################################################
# Catalog construction
################################################################################

def test_directory_from_environment(tmp_path: Path,
                                    monkeypatch: pytest.MonkeyPatch) -> None:
    write_ybsc_catalog(tmp_path, [VEGA])
    monkeypatch.setenv('YBSC_PATH', str(tmp_path))

    cat = YBSCStarCatalog()

    assert len(cat._stars) == 1
    assert cat.debug_level == 0


def test_missing_directory_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv('YBSC_PATH', raising=False)

    with pytest.raises(KeyError):
        YBSCStarCatalog()


def test_directory_argument_types(tmp_path: Path) -> None:
    write_ybsc_catalog(tmp_path, [VEGA])

    for dir_arg in (str(tmp_path), tmp_path, FCPath(tmp_path)):
        assert len(YBSCStarCatalog(dir_arg)._stars) == 1


def test_unknown_field_rejected() -> None:
    with pytest.raises(KeyError):
        ybsc_record(not_a_field='1')

    with pytest.raises(ValueError):
        ybsc_record(hr='123456')


################################################################################
# YBSCStar.__str__
################################################################################

def test_str_populated(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [VEGA])

    s = str(cat._stars[0])

    assert 'UNIQUE ID 7001' in s
    assert 'Name "3Alp Lyr"' in s
    assert 'Durch "BD+38 3238"' in s
    assert 'Draper 172167' in s
    assert 'SAO 67174' in s
    assert 'FK5 699' in s
    assert 'IR 1 Ref NASA' in s
    assert 'Multiple " "' in s
    assert 'Aitken 11510 AB' in s
    assert 'Variable "Alp Lyr"' in s
    assert 'SCLASS Code v' in s
    assert 'Galactic LON 67.44' in s
    assert 'LAT 19.24' in s
    assert 'B-V 0.0 | U-B -0.01 | R-I -0.03' in s
    assert 'Parallax DYN 0.1230000 arcsec' in s
    assert 'RadVel -14.0 km/s' in s
    assert 'RotVel (v sin i) 15.0 km/s' in s
    assert 'Double mag diff 10.40 Sep 62.80 arcsec' in s
    assert 'Components AB' in s
    assert '# 5' in s


def test_str_unfilled(tmp_path: Path) -> None:
    minimal = ybsc_record(hr='1', ir_ref='X', **ybsc_ra(0, 0, 1.1),
                          **ybsc_dec('+', 44, 40, 22),
                          glon='114.44', glat='-16.88', vmag='6.70',
                          sptype='A1Vn', pm_ra='-0.012', pm_dec='-0.018')
    cat = make_catalog(tmp_path, [minimal])
    star = cat._stars[0]
    star.galactic_longitude = None
    star.galactic_latitude = None

    s = str(star)

    assert 'IR 0 Ref N/A' in s
    assert 'Galactic LON N/A LAT N/A' in s
    assert 'Parallax TRIG N/A' in s
    assert 'RadVel N/A' in s
    assert 'RotVel (v sin i) N/A' in s
    assert 'Double mag diff N/A Sep N/A' in s
    assert '# N/A' in s
