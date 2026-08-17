################################################################################
# tests/test_spice.py
################################################################################

from __future__ import annotations

from pathlib import Path
from typing import Any

import cspyce
import numpy as np
import pytest
from catalog_data import SpiceTestStar, write_spice_catalog
from filecache import FCPath

from starcat import SpiceStar, SpiceStarCatalog, Star

# RA and DEC are in degrees
TEST_STARS = [SpiceTestStar(10., -10., 1., 'G2', 101),
              SpiceTestStar(20., 0., 5., 'A0*', 102),
              SpiceTestStar(30., 10., 9., 'M5', 103),
              SpiceTestStar(350., 20., 13., 'QQ', 104)]


@pytest.fixture(scope='session')
def catalog_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Write the star catalog shared by most of these tests."""

    root = tmp_path_factory.mktemp('spice')
    write_spice_catalog(root, 'testcat', TEST_STARS)

    return root


@pytest.fixture
def cat(catalog_dir: Path) -> SpiceStarCatalog:
    return SpiceStarCatalog('testcat', dir=catalog_dir)


def numbers(stars: list[Star]) -> list[int | None]:
    return [star.unique_number for star in stars]


################################################################################
# SpiceStar
################################################################################

def test_spice_star_defaults() -> None:
    star = SpiceStar()

    assert star.unique_number is None
    assert star.ra is None
    assert star.ra_sigma is None
    assert star.dec is None
    assert star.dec_sigma is None
    assert star.vmag is None
    assert star.spectral_class is None
    assert star.temperature is None


################################################################################
# Reading records
################################################################################

def test_parse_star(cat: SpiceStarCatalog) -> None:
    stars = list(cat.find_stars())
    assert len(stars) == 4

    star = stars[0]
    assert isinstance(star, SpiceStar)
    assert star.unique_number == 101
    assert star.ra == pytest.approx(np.radians(10.))
    assert star.dec == pytest.approx(np.radians(-10.))
    assert star.ra_sigma == pytest.approx(np.radians(0.001))
    assert star.dec_sigma == pytest.approx(np.radians(0.002))
    assert star.vmag == pytest.approx(1.)
    assert star.spectral_class == 'G2'
    assert star.temperature == 5780


def test_spectral_class_with_asterisk(cat: SpiceStarCatalog) -> None:
    star = next(s for s in cat.find_stars() if s.unique_number == 102)

    assert star.spectral_class == 'A0*'
    assert star.temperature == 10800


def test_unknown_spectral_class(cat: SpiceStarCatalog) -> None:
    star = next(s for s in cat.find_stars() if s.unique_number == 104)

    assert star.spectral_class == 'QQ'
    assert star.temperature is None


def test_full_result_false(cat: SpiceStarCatalog) -> None:
    stars = list(cat.find_stars(full_result=False))

    assert len(stars) == 4
    for star in stars:
        assert star.ra is not None
        assert star.vmag is not None
        assert star.temperature is None


def test_count_stars(cat: SpiceStarCatalog) -> None:
    assert cat.count_stars() == 4


################################################################################
# Searching
################################################################################

def test_find_stars_by_position(cat: SpiceStarCatalog) -> None:
    stars = list(cat.find_stars(ra_min=np.radians(15.), ra_max=np.radians(35.),
                                dec_min=np.radians(-5.),
                                dec_max=np.radians(15.)))

    assert numbers(stars) == [102, 103]


def test_find_stars_ra_wraparound(cat: SpiceStarCatalog) -> None:
    stars = list(cat.find_stars(ra_min=np.radians(340.),
                                ra_max=np.radians(15.)))

    assert set(numbers(stars)) == {101, 104}


def test_find_stars_by_vmag(cat: SpiceStarCatalog) -> None:
    assert numbers(list(cat.find_stars(vmag_min=5.))) == [102, 103, 104]
    assert numbers(list(cat.find_stars(vmag_max=9.))) == [101, 102, 103]
    assert numbers(list(cat.find_stars(vmag_min=5., vmag_max=9.))) == [102, 103]
    assert cat.count_stars(vmag_max=10.) == 3


def test_star_without_magnitude(cat: SpiceStarCatalog,
                                monkeypatch: pytest.MonkeyPatch) -> None:
    # SPICE always supplies a magnitude, but the reader tolerates one that is
    # missing by ignoring the magnitude limits
    def fake_stcg01(index: int) -> tuple[Any, ...]:
        return (0.1, 0.2, 1e-6, 2e-6, 42, 'G2', None)

    monkeypatch.setattr(cspyce, 'stcg01', fake_stcg01)

    stars = list(cat.find_stars(vmag_min=100., vmag_max=200.))

    assert len(stars) == 4
    assert stars[0].vmag is None
    assert stars[0].unique_number == 42


def test_debug_output(cat: SpiceStarCatalog,
                      capsys: pytest.CaptureFixture[str]) -> None:
    cat.debug_level = 1

    assert numbers(list(cat.find_stars(vmag_min=5., vmag_max=9.))) == [102, 103]

    out = capsys.readouterr().out
    assert out.count('SKIPPED VMAG') == 2
    assert 'SKIPPED VMAG 1.0' in out
    assert 'SKIPPED VMAG 13.0' in out
    assert 'UNIQUE ID 102' in out
    assert 'UNIQUE ID 103' in out
    assert '-' * 80 in out


################################################################################
# Catalog construction
################################################################################

def test_directory_argument_types(catalog_dir: Path) -> None:
    for dir_arg in (str(catalog_dir), catalog_dir, FCPath(catalog_dir)):
        cat = SpiceStarCatalog('testcat', dir=dir_arg)
        assert cat.count_stars() == 4
        assert cat.debug_level == 0


def test_spice_path_environment(tmp_path: Path,
                                monkeypatch: pytest.MonkeyPatch) -> None:
    write_spice_catalog(tmp_path / 'Stars', 'spicepathcat', TEST_STARS)
    monkeypatch.setenv('SPICE_PATH', str(tmp_path))
    monkeypatch.delenv('OOPS_RESOURCES', raising=False)

    cat = SpiceStarCatalog('spicepathcat')

    assert cat.count_stars() == 4


def test_oops_resources_environment(tmp_path: Path,
                                    monkeypatch: pytest.MonkeyPatch) -> None:
    write_spice_catalog(tmp_path / 'SPICE' / 'Stars', 'oopscat', TEST_STARS)
    monkeypatch.delenv('SPICE_PATH', raising=False)
    monkeypatch.setenv('OOPS_RESOURCES', str(tmp_path))

    cat = SpiceStarCatalog('oopscat')

    assert cat.count_stars() == 4


def test_no_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv('SPICE_PATH', raising=False)
    monkeypatch.delenv('OOPS_RESOURCES', raising=False)

    with pytest.raises(RuntimeError):
        SpiceStarCatalog('testcat')


def test_missing_catalog_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        SpiceStarCatalog('does_not_exist', dir=tmp_path)
