################################################################################
# tests/test_ucac4.py
################################################################################

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from filecache import FCPath
import numpy as np
import pytest

from starcat import Star, UCAC4Star, UCAC4StarCatalog
from starcat.starcatalog import HALFPI, MAS_TO_RAD, TWOPI, YEAR_TO_SEC
from starcat.ucac4 import UCAC4_RECORD_SIZE

from catalog_data import ucac4_ra, ucac4_record, ucac4_spd, write_ucac4_zone


# Zone 451 covers declinations 0.0 to 0.2 degrees
ZNUM = 451
DEC_MIN = np.radians(0.)
DEC_MAX = np.radians(0.2)
DEC_DEFAULT = 0.1
COS_DEC = np.cos(np.radians(DEC_DEFAULT))


def record(ra_deg: float = 10.,
           dec_deg: float = DEC_DEFAULT,
           **fields: int) -> bytes:
    """Return a packed UCAC4 record at the given RA and DEC in degrees."""

    return ucac4_record(ra=ucac4_ra(ra_deg), spd=ucac4_spd(dec_deg), **fields)


def make_catalog(root: Path,
                 records: Sequence[bytes],
                 znum: int = ZNUM) -> UCAC4StarCatalog:
    """Write a one-zone UCAC4 catalog and return a catalog object for it."""

    write_ucac4_zone(root, znum, records)

    return UCAC4StarCatalog(root)


def find(cat: UCAC4StarCatalog, **kwargs: Any) -> list[Star]:
    """Return the stars in the zone used by most of these tests."""

    return list(cat.find_stars(dec_min=DEC_MIN, dec_max=DEC_MAX, **kwargs))


################################################################################
# UCAC4Star
################################################################################

def test_ucac4_star_defaults() -> None:
    star = UCAC4Star()

    # Inherited from Star
    assert star.unique_number is None
    assert star.ra is None
    assert star.dec is None

    # UCAC4-specific
    assert star.vmag_model is None
    assert star.obj_type is None
    assert star.double_star_flag is None
    assert star.double_star_type is None
    assert star.galaxy_match is None
    assert star.extended_source is None
    assert star.num_img_total is None
    assert star.num_img_used is None
    assert star.num_cat_pm is None
    assert star.ra_mean_epoch is None
    assert star.dec_mean_epoch is None
    assert star.cat_match is None
    assert star.apass_mag_b is None
    assert star.apass_mag_v is None
    assert star.apass_mag_g is None
    assert star.apass_mag_r is None
    assert star.apass_mag_i is None
    assert star.apass_mag_b_sigma is None
    assert star.apass_mag_v_sigma is None
    assert star.apass_mag_g_sigma is None
    assert star.apass_mag_r_sigma is None
    assert star.apass_mag_i_sigma is None
    assert star.johnson_mag_b is None
    assert star.johnson_mag_v is None
    assert star.id_str is None
    assert star.id_str_ucac2 is None


def test_ucac4_star_str_empty() -> None:
    star = UCAC4Star()

    s = str(star)

    assert 'UNIQUE ID None' in s
    assert 'OBJTYPE None' in s
    assert 'APER VMAG None' in s
    assert 'DBL STAR FLAG=None TYPE=None' in s
    assert 'GALAXY NONE' in s
    assert 'EXT SOURCE NONE' in s
    assert 'APASS B None V None G None R None I None' in s
    assert 'JOHNSON B None V None' in s


def test_ucac4_star_str_populated(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path,
                       [record(cdf=23, leda=3, x2m=5,
                               apasm_g=10800, apasm_r=10900, apasm_i=11100,
                               apase_g=4, apase_r=5, apase_i=6)])

    star = find(cat, allow_galaxy=True)[0]
    s = str(star)

    assert 'OBJTYPE CLEAN' in s
    assert 'APER VMAG 10.000' in s
    assert 'DBL STAR FLAG=COMP2 TYPE=SECONDARY_PEAK' in s
    assert 'GALAXY Yes' in s
    assert 'EXT SOURCE Yes' in s
    assert ('APASS B 11.000 +/-  0.020 V 10.500 +/-  0.030 '
            'G 10.800 +/-  0.040 R 10.900 +/-  0.050 '
            'I 11.100 +/-  0.060') in s
    assert 'JOHNSON B 11.000 V 10.500' in s


def test_ucac4_star_str_no_galaxy(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record()])

    s = str(find(cat)[0])

    assert 'GALAXY No' in s
    assert 'EXT SOURCE No' in s


################################################################################
# Reading records
################################################################################

def test_parse_star(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., dec_deg=DEC_DEFAULT)])

    stars = find(cat)
    assert len(stars) == 1
    star = stars[0]

    assert isinstance(star, UCAC4Star)

    # Position
    assert star.ra == pytest.approx(np.radians(10.))
    assert star.dec == pytest.approx(np.radians(DEC_DEFAULT))
    assert star.rac_sigma == pytest.approx(20*MAS_TO_RAD)
    assert star.ra_sigma == pytest.approx(20*MAS_TO_RAD/COS_DEC)
    assert star.dec_sigma == pytest.approx(30*MAS_TO_RAD)

    # Magnitudes
    assert star.vmag_model == pytest.approx(10.)
    assert star.vmag == pytest.approx(10.5)
    assert star.vmag_sigma == pytest.approx(0.05)

    # Flags
    assert star.obj_type == UCAC4Star.UCAC4_OBJ_TYPE_CLEAN
    assert star.double_star_flag == UCAC4Star.UCAC4_DOUBLE_STAR_FLAG_SINGLE
    assert star.double_star_type == UCAC4Star.UCAC4_DOUBLE_STAR_TYPE_NONE
    assert star.galaxy_match == 0
    assert star.extended_source == 0

    # Images and epochs
    assert star.num_img_total == 6
    assert star.num_img_used == 5
    assert star.num_cat_pm == 2
    assert star.ra_mean_epoch == pytest.approx(1995.)
    assert star.dec_mean_epoch == pytest.approx(1996.)

    # Proper motion
    assert star.pm_rac == pytest.approx(10*MAS_TO_RAD*YEAR_TO_SEC)
    assert star.pm_ra == pytest.approx(10*MAS_TO_RAD*YEAR_TO_SEC/COS_DEC)
    assert star.pm_dec == pytest.approx(-5*MAS_TO_RAD*YEAR_TO_SEC)
    assert star.pm_rac_sigma == pytest.approx(2*MAS_TO_RAD*YEAR_TO_SEC)
    assert star.pm_ra_sigma == pytest.approx(2*MAS_TO_RAD*YEAR_TO_SEC/COS_DEC)
    assert star.pm_dec_sigma == pytest.approx(1*MAS_TO_RAD*YEAR_TO_SEC)

    # Photometry
    assert star.apass_mag_b == pytest.approx(11.)
    assert star.apass_mag_v == pytest.approx(10.5)
    assert star.apass_mag_g is None
    assert star.apass_mag_r is None
    assert star.apass_mag_i is None
    assert star.apass_mag_b_sigma == pytest.approx(0.02)
    assert star.apass_mag_v_sigma == pytest.approx(0.03)
    assert star.apass_mag_g_sigma == pytest.approx(0.99)
    assert star.apass_mag_r_sigma == pytest.approx(0.99)
    assert star.apass_mag_i_sigma == pytest.approx(0.99)
    assert star.johnson_mag_b == pytest.approx(11.)
    assert star.johnson_mag_v == pytest.approx(10.5)

    # Derived spectral class and temperature; B-V = 0.5 is closest to F7
    assert star.spectral_class == 'F7'
    assert star.temperature == 6313

    # Identification
    assert star.cat_match == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert star.unique_number == 1000000
    assert star.id_str == 'UCAC4-451-000001'
    assert star.id_str_ucac2 is None


def test_parse_negative_dec(tmp_path: Path) -> None:
    # Zone 1 covers declinations -90.0 to -89.8 degrees
    cat = make_catalog(tmp_path, [record(dec_deg=-89.9)], znum=1)

    stars = list(cat.find_stars(dec_min=np.radians(-90.),
                                dec_max=np.radians(-89.8)))
    assert len(stars) == 1
    assert stars[0].dec == pytest.approx(np.radians(-89.9))
    assert stars[0].id_str == 'UCAC4-001-000001'


def test_ucac2_id(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., zn2=0, rn2=0),
                                  record(ra_deg=20., zn2=123, rn2=456)])

    stars = find(cat)

    assert stars[0].id_str_ucac2 is None
    assert stars[1].id_str_ucac2 == 'UCAC2-123-000456'
    assert stars[0].id_str == 'UCAC4-451-000001'
    assert stars[1].id_str == 'UCAC4-451-000002'


def test_missing_magnitudes(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(magm=20000, maga=20000, sigmag=99)])

    star = find(cat)[0]

    assert star.vmag_model is None
    assert star.vmag is None
    assert star.vmag_sigma is None


def test_missing_photometry(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(apasm_b=20000, apasm_v=20000,
                                         apase_b=99, apase_v=99)])

    star = find(cat)[0]

    assert star.apass_mag_b is None
    assert star.apass_mag_v is None
    assert star.johnson_mag_b is None
    assert star.johnson_mag_v is None
    assert star.spectral_class is None
    assert star.temperature is None


def test_tycho_photometry_conversion(tmp_path: Path) -> None:
    # An APASS error of zero means the magnitudes are really Tycho Bt and Vt
    cat = make_catalog(tmp_path, [record(apasm_b=11000, apasm_v=10500,
                                         apase_b=0, apase_v=0)])

    star = find(cat)[0]

    assert star.apass_mag_b == pytest.approx(11.)
    assert star.apass_mag_v == pytest.approx(10.5)
    assert star.johnson_mag_v == pytest.approx(10.5 - 0.09*0.5)
    assert star.johnson_mag_b == pytest.approx(10.5 - 0.09*0.5 + 0.85*0.5)
    # B-V is now 0.425, which is closest to F5
    assert star.spectral_class == 'F5'
    assert star.temperature == 6540


def test_no_conversion_when_only_one_error_is_zero(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(apase_b=0, apase_v=3)])

    star = find(cat)[0]

    assert star.johnson_mag_b == pytest.approx(11.)
    assert star.johnson_mag_v == pytest.approx(10.5)


def test_galaxy_match_size(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(leda=3)])

    star = find(cat, allow_galaxy=True)[0]

    # log10 of the diameter in units of 0.1 arcmin, converted to degrees
    assert star.galaxy_match == pytest.approx(1000./10./60.)


################################################################################
# Proper motion special cases
################################################################################

def test_pm_too_large(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., pmrac=32767),
                                  record(ra_deg=20., pmdc=32767)])

    stars = find(cat, require_pm=False)

    for star in stars:
        assert star.pm_rac is None
        assert star.pm_ra is None
        assert star.pm_dec is None


@pytest.mark.parametrize('code,expected',
                         [(251, 27.5), (252, 32.5), (253, 37.5), (254, 45.)])
def test_pm_ra_sigma_codes(tmp_path: Path,
                           code: int,
                           expected: float) -> None:
    cat = make_catalog(tmp_path, [record(sigpmr=code-128)])

    star = find(cat)[0]

    assert star.pm_rac_sigma == pytest.approx(expected*MAS_TO_RAD*YEAR_TO_SEC)
    assert star.pm_ra_sigma == pytest.approx(expected*MAS_TO_RAD*YEAR_TO_SEC /
                                             COS_DEC)


@pytest.mark.parametrize('code,expected',
                         [(251, 27.5), (252, 32.5), (253, 37.5), (254, 45.)])
def test_pm_dec_sigma_codes(tmp_path: Path,
                            code: int,
                            expected: float) -> None:
    cat = make_catalog(tmp_path, [record(sigpmd=code-128)])

    star = find(cat)[0]

    assert star.pm_dec_sigma == pytest.approx(expected*MAS_TO_RAD)


def test_pm_no_data(tmp_path: Path) -> None:
    # A sigma of 255 with a proper motion of zero means "no data"
    cat = make_catalog(tmp_path, [record(pmrac=0, pmdc=0,
                                         sigpmr=255-128, sigpmd=255-128)])

    stars = find(cat, require_pm=False)
    assert len(stars) == 1
    star = stars[0]

    assert star.pm_rac is None
    assert star.pm_ra is None
    assert star.pm_rac_sigma is None
    assert star.pm_ra_sigma is None
    assert star.pm_dec is None
    assert star.pm_dec_sigma is None

    # Such stars are excluded by default
    assert find(cat) == []


def test_pm_no_sigma_but_nonzero_pm(tmp_path: Path) -> None:
    # A sigma of 255 with a non-zero proper motion keeps the proper motion
    cat = make_catalog(tmp_path, [record(sigpmr=255-128, sigpmd=255-128)])

    star = find(cat)[0]

    assert star.pm_rac == pytest.approx(10*MAS_TO_RAD*YEAR_TO_SEC)
    assert star.pm_dec == pytest.approx(-5*MAS_TO_RAD*YEAR_TO_SEC)
    assert star.pm_rac_sigma is None
    assert star.pm_ra_sigma is None
    assert star.pm_dec_sigma is None


################################################################################
# Selection criteria
################################################################################

def test_vmag_limits(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., maga=5000),
                                  record(ra_deg=20., maga=10000),
                                  record(ra_deg=30., maga=15000),
                                  record(ra_deg=40., maga=20000)])  # No data

    assert len(find(cat)) == 4
    assert [s.vmag for s in find(cat, vmag_min=8.)] == [10., 15.]
    assert [s.vmag for s in find(cat, vmag_max=12.)] == [5., 10.]
    assert [s.vmag for s in find(cat, vmag_min=8., vmag_max=12.)] == [10.]


@pytest.mark.parametrize('obj_type,is_clean',
                         [(0, True), (1, True), (2, False), (3, True),
                          (4, True), (5, True), (6, True), (7, True),
                          (8, False), (9, False)])
def test_require_clean(tmp_path: Path,
                       obj_type: int,
                       is_clean: bool) -> None:
    cat = make_catalog(tmp_path, [record(objt=obj_type)])

    assert len(find(cat)) == (1 if is_clean else 0)
    assert len(find(cat, require_clean=False)) == 1
    assert find(cat, require_clean=False)[0].obj_type == obj_type


def test_allow_double(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., cdf=0),
                                  record(ra_deg=20., cdf=23)])

    assert len(find(cat)) == 2
    assert len(find(cat, allow_double=False)) == 1

    double = find(cat)[1]
    assert double.double_star_flag == UCAC4Star.UCAC4_DOUBLE_STAR_FLAG_COMP2
    assert (double.double_star_type ==
            UCAC4Star.UCAC4_DOUBLE_STAR_TYPE_SECONDARY_PEAK)


def test_allow_galaxy(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10.),
                                  record(ra_deg=20., leda=2),
                                  record(ra_deg=30., x2m=7)])

    assert len(find(cat)) == 1
    assert len(find(cat, allow_galaxy=True)) == 3
    assert find(cat, allow_galaxy=True)[2].extended_source == 7


def test_require_pm(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10.),
                                  record(ra_deg=20., pmrac=0, pmdc=0,
                                         sigpmr=255-128, sigpmd=255-128)])

    assert len(find(cat)) == 1
    assert len(find(cat, require_pm=False)) == 2


def test_return_everything(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., objt=2),
                                  record(ra_deg=20., cdf=23),
                                  record(ra_deg=30., leda=2),
                                  record(ra_deg=40., pmrac=0, pmdc=0,
                                         sigpmr=255-128, sigpmd=255-128)])

    assert len(find(cat, allow_double=False)) == 0
    assert len(find(cat, return_everything=True)) == 4

    # return_everything overrides the individual options
    assert len(find(cat, return_everything=True, require_clean=True,
                    allow_double=False, allow_galaxy=False,
                    require_pm=True)) == 4


def test_ra_dec_limits(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., dec_deg=0.05),
                                  record(ra_deg=20., dec_deg=0.10),
                                  record(ra_deg=30., dec_deg=0.15)])

    # The limits are computed exactly the way the reader computes the stored
    # values so that the inclusive/exclusive behavior can be pinned down
    def ra_limit(ra_deg: float) -> float:
        return float(ucac4_ra(ra_deg) * MAS_TO_RAD)

    def dec_limit(dec_deg: float) -> float:
        return float(ucac4_spd(dec_deg) * MAS_TO_RAD - HALFPI)

    # RA is inclusive at the bottom and exclusive at the top
    stars = list(cat.find_stars(ra_min=ra_limit(20.), ra_max=ra_limit(30.),
                                dec_min=DEC_MIN, dec_max=DEC_MAX))
    assert [s.id_str for s in stars] == ['UCAC4-451-000002']

    # DEC is likewise inclusive at the bottom and exclusive at the top
    stars = list(cat.find_stars(dec_min=dec_limit(0.10),
                                dec_max=dec_limit(0.15)))
    assert [s.id_str for s in stars] == ['UCAC4-451-000002']


def test_full_result_false(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record()])

    star = find(cat, full_result=False)[0]

    # The selection criteria are still filled in
    assert star.ra is not None
    assert star.dec is not None
    assert star.vmag is not None
    assert star.obj_type is not None
    assert star.pm_ra is not None

    # ...but the expensive extras are not
    assert star.unique_number is None
    assert star.ra_sigma is None
    assert star.dec_sigma is None
    assert star.num_img_total is None
    assert star.apass_mag_b is None
    assert star.johnson_mag_b is None
    assert star.cat_match is None
    assert star.id_str is None
    assert star.spectral_class is None
    assert star.temperature is None


def test_count_stars(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10.),
                                  record(ra_deg=20.),
                                  record(ra_deg=30., objt=2)])

    assert cat.count_stars(dec_min=DEC_MIN, dec_max=DEC_MAX) == 2
    assert cat.count_stars(dec_min=DEC_MIN, dec_max=DEC_MAX,
                           require_clean=False) == 3


################################################################################
# Zone files
################################################################################

def test_zone_filename(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record()])

    assert cat._zone_filename(1) == FCPath(tmp_path) / 'u4b' / 'z001'
    assert cat._zone_filename(451) == FCPath(tmp_path) / 'u4b' / 'z451'
    assert cat._zone_filename(900) == FCPath(tmp_path) / 'u4b' / 'z900'


def test_multiple_zones(tmp_path: Path) -> None:
    # Zone 451 is DEC 0.0 to 0.2 and zone 452 is DEC 0.2 to 0.4
    write_ucac4_zone(tmp_path, 451, [record(ra_deg=10., dec_deg=0.1)])
    write_ucac4_zone(tmp_path, 452, [record(ra_deg=20., dec_deg=0.3)])
    cat = UCAC4StarCatalog(tmp_path)

    stars = list(cat.find_stars(dec_min=np.radians(0.), dec_max=np.radians(0.4)))

    assert [s.id_str for s in stars] == ['UCAC4-451-000001', 'UCAC4-452-000001']


def test_zone_boundary_opens_no_extra_zone(tmp_path: Path) -> None:
    # Only zone 451 exists here, so a search whose dec_max is exactly the top
    # of zone 451 must not try to open zone 452, which could contain no star
    # inside the search box anyway
    cat = make_catalog(tmp_path, [record(dec_deg=0.1)])

    stars = list(cat.find_stars(dec_min=np.radians(0.), dec_max=np.radians(0.2)))

    assert [s.id_str for s in stars] == ['UCAC4-451-000001']


def test_zone_900(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(dec_deg=89.9)], znum=900)

    stars = list(cat.find_stars(dec_min=np.radians(89.8),
                                dec_max=np.radians(90.)))

    assert [s.id_str for s in stars] == ['UCAC4-900-000001']


def test_directory_from_environment(tmp_path: Path,
                                    monkeypatch: pytest.MonkeyPatch) -> None:
    write_ucac4_zone(tmp_path, ZNUM, [record()])
    monkeypatch.setenv('UCAC4_PATH', str(tmp_path))

    cat = UCAC4StarCatalog()

    assert len(find(cat)) == 1


def test_missing_directory_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv('UCAC4_PATH', raising=False)

    with pytest.raises(KeyError):
        UCAC4StarCatalog()


################################################################################
# RA optimization
################################################################################

def ra_search_records() -> list[bytes]:
    """Return records at 10, 20, 20, 30, and 40 degrees of RA."""

    return [record(ra_deg=10.), record(ra_deg=20.), record(ra_deg=20.),
            record(ra_deg=30.), record(ra_deg=40.)]


@pytest.mark.parametrize('ra_min_deg,expected',
                         [(0., [1, 2, 3, 4, 5]),
                          (5., [1, 2, 3, 4, 5]),
                          (10., [1, 2, 3, 4, 5]),
                          (15., [2, 3, 4, 5]),
                          (20., [2, 3, 4, 5]),     # Duplicated RA
                          (25., [4, 5]),
                          (40., [5]),
                          (45., [])])
def test_optimize_ra(tmp_path: Path,
                     ra_min_deg: float,
                     expected: list[int]) -> None:
    cat = make_catalog(tmp_path, ra_search_records())

    kwargs: dict[str, Any] = {'dec_min': DEC_MIN, 'dec_max': DEC_MAX,
                              'ra_min': np.radians(ra_min_deg)}
    optimized = list(cat.find_stars(**kwargs))
    unoptimized = list(cat.find_stars(optimize_ra=False, **kwargs))

    # The record numbers in the IDs must not depend on the search method
    ids = [f'UCAC4-451-{num:06d}' for num in expected]
    assert [s.id_str for s in optimized] == ids
    assert [s.id_str for s in unoptimized] == ids


def test_find_starting_ra(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, ra_search_records())

    with cat._zone_filename(ZNUM).open(mode='rb') as fp:
        # Not worth searching
        assert cat._find_starting_ra(fp, 0.) == 0
        assert cat._find_starting_ra(fp, -1.) == 0

        # Before, at, between, and after the records
        assert cat._find_starting_ra(fp, np.radians(5.)) == 0
        assert cat._find_starting_ra(fp, np.radians(10.)) == 0
        assert cat._find_starting_ra(fp, np.radians(15.)) == 1
        assert cat._find_starting_ra(fp, np.radians(20.)) == 1
        assert cat._find_starting_ra(fp, np.radians(25.)) == 3
        assert cat._find_starting_ra(fp, np.radians(40.)) == 4
        assert cat._find_starting_ra(fp, np.radians(45.)) == 5

        # The file is left positioned at the record that was found
        cat._find_starting_ra(fp, np.radians(25.))
        assert fp.tell() == 3 * UCAC4_RECORD_SIZE


def test_find_starting_ra_exact(tmp_path: Path) -> None:
    # An RA that is exactly representable in the file's units
    ra_mas = ucac4_ra(20.)
    cat = make_catalog(tmp_path, ra_search_records())

    with cat._zone_filename(ZNUM).open(mode='rb') as fp:
        assert cat._find_starting_ra(fp, ra_mas*MAS_TO_RAD) == 1


def test_truncated_zone_file(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record()[:-1]])

    with pytest.raises(AssertionError):
        find(cat, ra_min=np.radians(5.))


################################################################################
# Debugging output
################################################################################

def test_debug_level_1(tmp_path: Path,
                       capsys: pytest.CaptureFixture[str]) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., rnm=1, objt=2),
                                  record(ra_deg=20., rnm=2, cdf=23),
                                  record(ra_deg=30., rnm=3, leda=2),
                                  record(ra_deg=40., rnm=4, pmrac=0, pmdc=0,
                                         sigpmr=255-128, sigpmd=255-128),
                                  record(ra_deg=50., rnm=5)])
    cat.debug_level = 1

    stars = find(cat, allow_double=False)
    assert len(stars) == 1

    out = capsys.readouterr().out
    assert 'ID 1 SKIPPED NOT CLEAN 2' in out
    assert 'ID 2 SKIPPED DOUBLE 23' in out
    assert 'ID 3 SKIPPED GALAXY/EXTENDED 2 0' in out
    assert 'ID 4 SKIPPED NO PM' in out
    assert 'ID 5 OK!' in out
    assert 'UNIQUE ID 5' in out


def test_debug_level_1_count(tmp_path: Path,
                             capsys: pytest.CaptureFixture[str]) -> None:
    cat = make_catalog(tmp_path, [record(rnm=42)])
    cat.debug_level = 1

    assert cat.count_stars(dec_min=DEC_MIN, dec_max=DEC_MAX) == 1

    out = capsys.readouterr().out
    assert 'ID 42 OK!' in out
    # The star itself is not printed when the result is not complete
    assert 'UNIQUE ID' not in out


def test_debug_level_2(tmp_path: Path,
                       capsys: pytest.CaptureFixture[str]) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=10., rnm=1, dec_deg=0.19),
                                  record(ra_deg=20., rnm=2, maga=5000),
                                  record(ra_deg=30., rnm=3, maga=15000),
                                  record(ra_deg=40., rnm=4)])
    cat.debug_level = 2

    stars = list(cat.find_stars(ra_min=np.radians(15.), ra_max=np.radians(35.),
                                dec_min=DEC_MIN, dec_max=np.radians(0.15),
                                vmag_min=8., vmag_max=12.))
    assert stars == []

    out = capsys.readouterr().out
    assert 'ID 2 SKIPPED MODEL MAG' in out
    assert 'ID 3 SKIPPED MODEL MAG' in out
    assert 'ID 4 SKIPPED RA AND REST OF FILE' in out
    # The first record is before ra_min once the search is not optimized
    stars = list(cat.find_stars(ra_min=np.radians(15.), ra_max=np.radians(35.),
                                dec_min=DEC_MIN, dec_max=np.radians(0.15),
                                optimize_ra=False))
    assert [s.id_str for s in stars] == ['UCAC4-451-000002',
                                         'UCAC4-451-000003']
    assert 'ID 1 SKIPPED RA/DEC' in capsys.readouterr().out


################################################################################
# Wrap-around searches
################################################################################

def test_ra_wraparound(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=1.), record(ra_deg=180.),
                                  record(ra_deg=359.)])

    stars = list(cat.find_stars(ra_min=np.radians(350.), ra_max=np.radians(10.),
                                dec_min=DEC_MIN, dec_max=DEC_MAX))

    assert [s.id_str for s in stars] == ['UCAC4-451-000001',
                                         'UCAC4-451-000003']


def test_full_sky_ra(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record(ra_deg=359.9999)])

    stars = list(cat.find_stars(ra_min=0., ra_max=TWOPI,
                                dec_min=DEC_MIN, dec_max=DEC_MAX))

    assert len(stars) == 1


def test_optional_zone_argument(tmp_path: Path) -> None:
    write_ucac4_zone(tmp_path, ZNUM, [record()])

    # A string, a Path, and an FCPath must all work
    for dir_arg in (str(tmp_path), tmp_path, FCPath(tmp_path)):
        cat = UCAC4StarCatalog(dir_arg)
        assert len(find(cat)) == 1


def test_debug_level_default(tmp_path: Path) -> None:
    cat = make_catalog(tmp_path, [record()])

    assert cat.debug_level == 0


def test_unknown_field_rejected() -> None:
    with pytest.raises(KeyError):
        ucac4_record(not_a_field=1)
