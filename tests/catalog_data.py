################################################################################
# tests/catalog_data.py
#
# Build tiny, synthetic star catalogs in the on-disk formats used by the real
# UCAC4, YBSC, and SPICE catalogs.
#
# The real catalogs are multi-gigabyte data sets hosted outside of this
# repository. Rather than requiring them, the tests write small catalogs
# containing only the records they need into a temporary directory. The format
# descriptions below were transcribed from the published catalog documentation
# and are intentionally independent of the ones in the ``starcat`` package so
# that the tests exercise the readers instead of mirroring them.
################################################################################

from __future__ import annotations

import struct
from collections.abc import Iterable, Sequence
from pathlib import Path

import cspyce

################################################################################
# UCAC4
#
# Each zone file is a headerless sequence of 78-byte little-endian binary
# records sorted by increasing RA. See readme_u4.txt in the UCAC4 distribution.
################################################################################

UCAC4_FMT = '<iihhbbbbbbbbhhhhbbihhhbbbbbbhhhhhbbbbbbibbihi'
UCAC4_RECORD_SIZE = 78

# The fields of a UCAC4 record, in the order they appear in UCAC4_FMT.
UCAC4_FIELDS: tuple[str, ...] = (
    'ra',        # I*4 mas       Right ascension at J2000
    'spd',       # I*4 mas       South pole distance at J2000
    'magm',      # I*2 millimag  Fit model magnitude
    'maga',      # I*2 millimag  Aperture magnitude
    'sigmag',    # I*1 1/100 mag Error of UCAC magnitude
    'objt',      # I*1           Object type
    'cdf',       # I*1           Combined double star flag
    'sigra',     # I*1 mas       s.e. at central epoch in RA (*cos DEC); +128
    'sigdc',     # I*1 mas       s.e. at central epoch in DEC; +128
    'na1',       # I*1           Total # of CCD images
    'nu1',       # I*1           # of CCD images used
    'cu1',       # I*1           # catalogs used for proper motions
    'cepra',     # I*2 0.01 yr   Central epoch for mean RA, minus 1900
    'cepdc',     # I*2 0.01 yr   Central epoch for mean DEC, minus 1900
    'pmrac',     # I*2 0.1 mas/yr Proper motion in RA*cos(DEC)
    'pmdc',      # I*2 0.1 mas/yr Proper motion in DEC
    'sigpmr',    # I*1 0.1 mas/yr s.e. of pmRA*cos(DEC); +128
    'sigpmd',    # I*1 0.1 mas/yr s.e. of pmDEC; +128
    'pts_key',   # I*4           2MASS unique star identifier
    'j_m',       # I*2 millimag  2MASS J magnitude
    'h_m',       # I*2 millimag  2MASS H magnitude
    'k_m',       # I*2 millimag  2MASS K_s magnitude
    'icqflg1',   # I*1           2MASS quality flag for J
    'icqflg2',   # I*1           2MASS quality flag for H
    'icqflg3',   # I*1           2MASS quality flag for K_s
    'e2mpho1',   # I*1 1/100 mag Error 2MASS J magnitude
    'e2mpho2',   # I*1 1/100 mag Error 2MASS H magnitude
    'e2mpho3',   # I*1 1/100 mag Error 2MASS K_s magnitude
    'apasm_b',   # I*2 millimag  B magnitude from APASS
    'apasm_v',   # I*2 millimag  V magnitude from APASS
    'apasm_g',   # I*2 millimag  g magnitude from APASS
    'apasm_r',   # I*2 millimag  r magnitude from APASS
    'apasm_i',   # I*2 millimag  i magnitude from APASS
    'apase_b',   # I*1 1/100 mag Error of B magnitude from APASS
    'apase_v',   # I*1 1/100 mag Error of V magnitude from APASS
    'apase_g',   # I*1 1/100 mag Error of g magnitude from APASS
    'apase_r',   # I*1 1/100 mag Error of r magnitude from APASS
    'apase_i',   # I*1 1/100 mag Error of i magnitude from APASS
    'gcflg',     # I*1           Yale SPM g-flag*10 c-flag
    'icf',       # I*4           Catalog match flags
    'leda',      # I*1           LEDA galaxy match flag
    'x2m',       # I*1           2MASS extended source flag
    'rnm',       # I*4           Unique star identification number
    'zn2',       # I*2           Zone number of UCAC2 (0 = no match)
    'rn2',       # I*4           Running record number along UCAC2 zone
)

assert struct.calcsize(UCAC4_FMT) == UCAC4_RECORD_SIZE
assert len(UCAC4_FIELDS) == len(struct.unpack(UCAC4_FMT, b'\x00' * UCAC4_RECORD_SIZE))

# A "boring" star: clean, single, not a galaxy, with proper motion and with
# APASS B and V photometry (B-V = 0.5, which maps to spectral class F7).
UCAC4_DEFAULTS: dict[str, int] = {
    'ra': 0,
    'spd': 0,
    'magm': 10000,
    'maga': 10500,
    'sigmag': 5,
    'objt': 0,
    'cdf': 0,
    'sigra': -108,    # 20 mas
    'sigdc': -98,     # 30 mas
    'na1': 6,
    'nu1': 5,
    'cu1': 2,
    'cepra': 9500,    # 1995.0
    'cepdc': 9600,    # 1996.0
    'pmrac': 100,     # 10 mas/yr
    'pmdc': -50,      # -5 mas/yr
    'sigpmr': -108,   # 2 mas/yr
    'sigpmd': -118,   # 1 mas/yr
    'pts_key': 1234567,
    'j_m': 9500,
    'h_m': 9300,
    'k_m': 9200,
    'icqflg1': 5,
    'icqflg2': 5,
    'icqflg3': 5,
    'e2mpho1': 2,
    'e2mpho2': 3,
    'e2mpho3': 3,
    'apasm_b': 11000,
    'apasm_v': 10500,
    'apasm_g': 20000,  # No data
    'apasm_r': 20000,  # No data
    'apasm_i': 20000,  # No data
    'apase_b': 2,
    'apase_v': 3,
    'apase_g': 99,     # No data
    'apase_r': 99,     # No data
    'apase_i': 99,     # No data
    'gcflg': 0,
    'icf': 123456789,
    'leda': 0,
    'x2m': 0,
    'rnm': 1000000,
    'zn2': 0,
    'rn2': 0,
}


def ucac4_ra(ra_deg: float) -> int:
    """Convert an RA in degrees to the UCAC4 integer field value (mas)."""

    return round(ra_deg * 3600000.)


def ucac4_spd(dec_deg: float) -> int:
    """Convert a DEC in degrees to the UCAC4 south pole distance field (mas)."""

    return round((dec_deg+90.) * 3600000.)


def ucac4_record(**fields: int) -> bytes:
    """Return one packed UCAC4 record, overriding the given fields.

    Any field not given takes its value from UCAC4_DEFAULTS.
    """

    values = dict(UCAC4_DEFAULTS)
    for name, value in fields.items():
        if name not in values:
            raise KeyError(f'Unknown UCAC4 field "{name}"')
        values[name] = value

    return struct.pack(UCAC4_FMT, *[values[name] for name in UCAC4_FIELDS])


def write_ucac4_zone(root: Path,
                     znum: int,
                     records: Iterable[bytes]) -> Path:
    """Write a UCAC4 zone file containing the given records.

    Parameters:
        root: The catalog root directory; the zone file is placed in ``u4b``.
        znum: The zone number (1 to 900).
        records: The already-packed records, which must be sorted by RA.

    Returns:
        The path of the zone file.
    """

    zone_dir = root / 'u4b'
    zone_dir.mkdir(parents=True, exist_ok=True)
    path = zone_dir / f'z{znum:03d}'
    path.write_bytes(b''.join(records))

    return path


################################################################################
# YBSC
#
# The catalog is a text file of fixed-width 197-character records, one star per
# line. See the ReadMe accompanying the Bright Star Catalogue, 5th Revised Ed.
################################################################################

YBSC_RECORD_SIZE = 197

# Field name -> (0-based start column, width)
YBSC_FIELDS: dict[str, tuple[int, int]] = {
    'hr': (0, 4),                   # Harvard Revised Number
    'name': (4, 10),                # Bayer and/or Flamsteed name
    'dm': (14, 11),                 # Durchmusterung identification
    'hd': (25, 6),                  # Henry Draper Catalog number
    'sao': (31, 6),                 # SAO Catalog number
    'fk5': (37, 4),                 # FK5 star number
    'ir_flag': (41, 1),             # I if infrared source
    'ir_ref': (42, 1),              # Coded reference for infrared source
    'multiple': (43, 1),            # Double or multiple-star code
    'ads': (44, 5),                 # Aitken's Double Star Catalog designation
    'ads_comp': (49, 2),            # ADS number components
    'var_id': (51, 9),              # Variable star identification
    'ra_dec_1900': (60, 15),        # RA and DEC at B1900 (unused by starcat)
    'ra_h': (75, 2),                # Hours RA, J2000
    'ra_m': (77, 2),                # Minutes RA, J2000
    'ra_s': (79, 4),                # Seconds RA, J2000
    'dec_sign_d': (83, 3),          # Sign and degrees DEC, J2000
    'dec_m': (86, 2),               # Arcminutes DEC, J2000
    'dec_s': (88, 2),               # Arcseconds DEC, J2000
    'glon': (90, 6),                # Galactic longitude
    'glat': (96, 6),                # Galactic latitude
    'vmag': (102, 5),               # Visual magnitude
    'vmag_code': (107, 1),          # Visual magnitude code
    'vmag_uncertainty': (108, 1),   # Uncertainty flag on V
    'b_v': (109, 5),                # B-V color
    'b_v_uncertainty': (114, 1),    # Uncertainty flag on B-V
    'u_b': (115, 5),                # U-B color
    'u_b_uncertainty': (120, 1),    # Uncertainty flag on U-B
    'r_i': (121, 5),                # R-I color
    'r_i_code': (126, 1),           # Code for R-I system
    'sptype': (127, 20),            # Spectral type
    'sptype_code': (147, 1),        # Spectral type code
    'pm_ra': (148, 6),              # Annual proper motion in RA (arcsec/yr)
    'pm_dec': (154, 6),             # Annual proper motion in DEC (arcsec/yr)
    'parallax_type': (160, 1),      # D for a dynamical parallax
    'parallax': (161, 5),           # Parallax (arcsec)
    'rad_vel': (166, 4),            # Heliocentric radial velocity (km/s)
    'rad_vel_comments': (170, 4),   # Radial velocity comments
    'rot_vel_limit': (174, 2),      # Rotational velocity limit characters
    'rot_vel': (176, 3),            # Rotational velocity v sin i (km/s)
    'rot_vel_uncertainty': (179, 1),  # Uncertainty/variability flag on RotVel
    'dmag': (180, 4),               # Magnitude difference of double
    'sep': (184, 6),                # Separation of components (arcsec)
    'mult_id': (190, 4),            # Identification of components in Dmag
    'mult_cnt': (194, 2),           # Number of components in a multiple
    'note_flag': (196, 1),          # * if there is a note
}


def ybsc_record(**fields: str) -> str:
    """Return one 197-character YBSC record containing the given fields.

    Every field not given is left blank. Values are right-justified within
    their columns, which is how the numeric fields of the real catalog are
    formatted; the reader strips whitespace in either case.
    """

    record = [' '] * YBSC_RECORD_SIZE
    for name, value in fields.items():
        if name not in YBSC_FIELDS:
            raise KeyError(f'Unknown YBSC field "{name}"')
        start, width = YBSC_FIELDS[name]
        if len(value) > width:
            raise ValueError(f'YBSC field "{name}" is limited to {width} '
                             f'characters but got "{value}"')
        record[start:start+width] = list(value.rjust(width))

    return ''.join(record)


def ybsc_ra(hours: int,
            minutes: int,
            seconds: float) -> dict[str, str]:
    """Return the YBSC RA fields for the given J2000 hours/minutes/seconds."""

    return {'ra_h': f'{hours:02d}',
            'ra_m': f'{minutes:02d}',
            'ra_s': f'{seconds:04.1f}'}


def ybsc_dec(sign: str,
             degrees: int,
             minutes: int,
             seconds: int) -> dict[str, str]:
    """Return the YBSC DEC fields for the given J2000 sign/deg/min/sec."""

    return {'dec_sign_d': f'{sign}{degrees:02d}',
            'dec_m': f'{minutes:02d}',
            'dec_s': f'{seconds:02d}'}


def write_ybsc_catalog(root: Path,
                       records: Iterable[str]) -> Path:
    """Write a YBSC ``catalog`` file containing the given records.

    Parameters:
        root: The catalog root directory.
        records: The already-formatted fixed-width records.

    Returns:
        The path of the catalog file.
    """

    root.mkdir(parents=True, exist_ok=True)
    path = root / 'catalog'
    path.write_text(''.join(f'{record}\n' for record in records))

    return path


################################################################################
# SPICE
#
# A SPICE type 1 star catalog is an EK (Event Kernel) containing a single table
# with seven specific columns, none of which may be declared NULLS_OK.
################################################################################

SPICE_COLUMNS = ('RA', 'DEC', 'RA_SIGMA', 'DEC_SIGMA', 'CATALOG_NUMBER',
                 'SPECTRAL_TYPE', 'VISUAL_MAGNITUDE')

SPICE_DECLS = ('DATATYPE = DOUBLE PRECISION, INDEXED = TRUE',
               'DATATYPE = DOUBLE PRECISION, INDEXED = TRUE',
               'DATATYPE = DOUBLE PRECISION, INDEXED = TRUE',
               'DATATYPE = DOUBLE PRECISION, INDEXED = TRUE',
               'DATATYPE = INTEGER, INDEXED = TRUE',
               'DATATYPE = CHARACTER*(4), INDEXED = TRUE',
               'DATATYPE = DOUBLE PRECISION, INDEXED = TRUE')


class SpiceTestStar:
    """One star to be written into a synthetic SPICE type 1 star catalog."""

    def __init__(self,
                 ra: float,
                 dec: float,
                 vmag: float,
                 spectral_class: str,
                 catalog_number: int,
                 ra_sigma: float = 0.001,
                 dec_sigma: float = 0.002) -> None:
        """Create a star; RA, DEC, and the uncertainties are in degrees."""

        self.ra = ra
        self.dec = dec
        self.vmag = vmag
        self.spectral_class = spectral_class
        self.catalog_number = catalog_number
        self.ra_sigma = ra_sigma
        self.dec_sigma = dec_sigma


def write_spice_catalog(root: Path,
                        name: str,
                        stars: Sequence[SpiceTestStar],
                        table_name: str | None = None) -> Path:
    """Write a SPICE type 1 star catalog containing the given stars.

    Parameters:
        root: The directory to write the catalog into.
        name: The catalog name; the file is called ``name.bdb``.
        stars: The stars to include.
        table_name: The name of the EK table; defaults to the catalog name in
            upper case. Each catalog loaded during a single test session needs
            a unique table name because SPICE looks stars up by table.

    Returns:
        The path of the catalog file.
    """

    root.mkdir(parents=True, exist_ok=True)
    path = root / f'{name}.bdb'
    path.unlink(missing_ok=True)

    if table_name is None:
        table_name = name.upper()

    nstars = len(stars)
    entszs = [1] * nstars
    nlflgs = [False] * nstars

    handle = cspyce.ekopn(str(path), name, 0)
    segno, rcptrs = cspyce.ekifld(handle, table_name, nstars,
                                  list(SPICE_COLUMNS), list(SPICE_DECLS))

    def add_d(column: str, values: list[float]) -> None:
        cspyce.ekacld(handle, segno, column, values, entszs, nlflgs, rcptrs)

    add_d('RA', [s.ra for s in stars])
    add_d('DEC', [s.dec for s in stars])
    add_d('RA_SIGMA', [s.ra_sigma for s in stars])
    add_d('DEC_SIGMA', [s.dec_sigma for s in stars])
    cspyce.ekacli(handle, segno, 'CATALOG_NUMBER',
                  [s.catalog_number for s in stars], entszs, nlflgs, rcptrs)
    cspyce.ekaclc(handle, segno, 'SPECTRAL_TYPE',
                  [s.spectral_class for s in stars], entszs, nlflgs, rcptrs)
    add_d('VISUAL_MAGNITUDE', [s.vmag for s in stars])

    cspyce.ekffld(handle, segno, rcptrs)
    cspyce.ekcls(handle)

    return path
