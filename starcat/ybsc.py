################################################################################
# starcat/ybsc.py
################################################################################

# The Bright Star Catalogue,  5th Revised Ed. (Preliminary Version)
#      Hoffleit D., Warren Jr W.H.
#     <Astronomical Data Center, NSSDC/ADC (1991)>
#     =1964BS....C......0H
# From ftp://cdsarc.u-strasbg.fr/cats/V/50/1
# http://tdc-www.harvard.edu/catalogs/bsc5.html

from __future__ import annotations

import numpy as np
import os
from pathlib import Path
from typing import Any, Iterator, Optional, cast

from filecache import FCPath

from .starcatalog import (AS_TO_RAD,
                          YEAR_TO_SEC,
                          Star,
                          StarCatalog
                          )


class YBSCStar(Star):
    """A holder for star attributes.

    This class includes attributes unique to the YBSC catalog.
    """

    YBSC_IR_NASA = 0
    YBSC_IR_ENGLES = 1
    YBSC_IR_UNCERTAIN = 2

    YBSC_IR_STRINGS = ['NASA', 'ENGLES', 'UNCERTAIN']

    YBSC_MULTIPLE_NONE = ' '
    YBSC_MULTIPLE_ASTROMETRIC = 'A'
    YBSC_MULTIPLE_DUPLICITY_OCCULTATION = 'D'
    YBSC_MULTIPLE_INNES = 'I'
    YBSC_MULTIPLE_ROSSITER = 'R'
    YBSC_MULTIPLE_DUPLICITY_SPECKLE = 'S'
    YBSC_MULTIPLE_WORLEY = 'W'

    YBSC_VMAG_UNCERTAINTY_V = ' '
    YBSC_VMAG_UNCERTAINTY_HR_REDUCED = 'R'
    YBSC_VMAG_UNCERTAINTY_HR = 'H'

    def __init__(self) -> None:
        # Initialize the standard fields
        super().__init__()

        # Initialize the YBSC-specific fields
        self.name: Optional[str] = None
        """Bayer and/or Flamsteed name"""

        self.durchmusterung_id: Optional[str] = None
        """Durchmusterung identification"""

        self.draper_number: Optional[int] = None
        """Henry Draper Catalog number (out of 225300)"""

        self.sao_number: Optional[int] = None
        """SAO Catalog number (out of 258997)"""

        self.fk5_number: Optional[int] = None
        """FK5 star number"""

        self.ir_source: Optional[bool] = None
        """True if infrared source"""

        self.ir_source_ref: Optional[int] = None
        """Infrared source:
           NASA, ENGLES, or UNCERTAIN"""

        self.multiple_star_code: Optional[str] = None
        """Double or multiple star code:
           'A' = Astrometric binary;
           'D' = Duplicity discovered by occultation;
           'I' = Innes, Southern Double Star Catalogue (1927);
           'R' = Rossiter, Michigan Publ. 9, 1955;
           'S' = Duplicity discovered by speckle interferometry;
           'W' = Worley (1978) update of the IDS"""

        self.aitken_designation: Optional[str] = None
        """Aitken's Double Star Catalog (ADS) designation"""

        self.ads_components: Optional[str] = None
        """ADS number components"""

        self.variable_star_id: Optional[str] = None
        """Variable star identification"""

        self.galactic_longitude: Optional[float] = None
        """Galactic longitude (radians)"""

        self.galactic_latitude: Optional[float] = None
        """Galactic latitude (radians)"""

        self.vmag_code: Optional[str] = None
        """Visual magnitude code:
           ' ' = V on UBV Johnson system;
           'R' = HR magnitudes reduced to the UBV system;
           'H' = original HR magnitude"""

        self.vmag_uncertainty_flag: Optional[str] = None
        """Uncertainty flag on visual magnitude"""

        self.b_v: Optional[float] = None
        """B-V color in the UBV system"""

        self.b_v_uncertainty_flag: Optional[str] = None
        """Uncertainty flag on B-V color"""

        self.u_b: Optional[float] = None
        """U-B color in the UBV system"""

        self.u_b_uncertainty_flag: Optional[str] = None
        """Uncertainty flag on U-B color"""

        self.r_i: Optional[float] = None
        """R-I color in the system indicated by r_i_code"""

        self.r_i_code: Optional[str] = None
        """Code for R-I system:
           'C' = Cousin;
           'E' = 'Eggen';
           ':' = Unknown;
           '?' = Unknown;
           'D' = Unknown"""

        self.spectral_class_code: Optional[str] = None
        """Spectral class code:
           'e', 'v', or 't'"""

        self.parallax_type: Optional[str] = None
        """Parallax type:
           'D' = Dyanmical, otherwise Trigonometric"""

        self.parallax: Optional[float] = None
        """Parallax (arcsec); see parallax_type for measurement type"""

        self.radial_velocity: Optional[float] = None
        """Radial velocity (km/s)"""

        self.radial_velocity_comments: Optional[str] = None
        """Radial velocity comments (multiple possible):
           'V' = Variable radial velocity;
           'V?' = Suspected variable radial velocity;
           'SB', 'SB1', 'SB2', 'SB3' = Spectroscopic binaries,
           single/double/triple-lined spectra;
           'O' = Orbital data available"""

        self.rotational_velocity_limit: Optional[str] = None
        """Rotational velocity limit:
           '<', '=', or '>'"""

        self.rotational_velocity: Optional[float] = None
        """Rotational velocity [v sin i] (km/s)"""

        self.rotational_velocity_uncertainty_flag: Optional[str] = None
        """Rotational velocity uncertainty and variability flag:
           ' ', ':', or 'v'"""

        self.double_mag_diff: Optional[float] = None
        """Magnitude difference of double, or brightest multiple"""

        self.double_mag_sep: Optional[float] = None
        """Separation of components in double_mag if occultation binary (radians)"""

        self.double_mag_components: Optional[str] = None
        """Indentification of components in double_mag"""

        self.multiple_num_components: Optional[int] = None
        """Number of components assigned to a multiple"""

    def __str__(self) -> str:
        ret = Star.__str__(self) + '\n'

        ret += f'Name "{self.name}"'
        ret += f' | Durch "{self.durchmusterung_id}"'
        ret += f' | Draper {self.draper_number}'
        ret += f' | SAO {self.sao_number}'
        ret += f' | FK5 {self.fk5_number}'
        ret += '\n'

        if self.ir_source_ref is not None:
            ir_ref = YBSCStar.YBSC_IR_STRINGS[self.ir_source_ref]
        else:
            ir_ref = 'N/A'
        ret += f'IR {self.ir_source:d} Ref {ir_ref}'
        ret += f' | Multiple "{self.multiple_star_code}"'
        ret += f' | Aitken {self.aitken_designation} '
        ret += str(self.ads_components)
        ret += f' | Variable "{self.variable_star_id}"'
        ret += '\n'

        ret += f'SCLASS Code {self.spectral_class_code}'

        ret += ' | Galactic LON '
        if self.galactic_longitude is None:
            ret += 'N/A'
        else:
            ret += str(np.degrees(self.galactic_longitude))
        ret += ' LAT '
        if self.galactic_latitude is None:
            ret += 'N/A'
        else:
            ret += str(np.degrees(self.galactic_latitude))
        ret += '\n'

        ret += f'B-V {self.b_v}'
        ret += f' | U-B {self.u_b}'
        ret += f' | R-I {self.r_i}'
        ret += '\n'

        ret += 'Parallax '
        if self.parallax_type == 'D':
            ret += 'DYN'
        else:
            ret += 'TRIG'

        if self.parallax is None:
            ret += ' N/A'
        else:
            ret += f' {self.parallax:.7f} arcsec'
        if self.radial_velocity is None:
            ret += ' | RadVel N/A'
        else:
            ret += f' | RadVel {self.radial_velocity} km/s'
        ret += f' {self.radial_velocity_comments}'
        ret += f' {self.rotational_velocity_limit}'
        if self.rotational_velocity is None:
            ret += ' | RotVel (v sin i) N/A'
        else:
            ret += f' | RotVel (v sin i) {self.rotational_velocity} km/s'
        ret += f' {self.rotational_velocity_uncertainty_flag}'
        ret += '\n'

        if self.double_mag_diff is None:
            ret += 'Double mag diff N/A'
        else:
            ret += f'Double mag diff {self.double_mag_diff:.2f}'
        if self.double_mag_sep is None:
            ret += ' Sep N/A'
        else:
            ret += f' Sep {self.double_mag_sep / AS_TO_RAD:.2f} arcsec'
        ret += f' Components {self.double_mag_components}'
        if self.multiple_num_components is None:
            ret += ' # N/A'
        else:
            ret += f' # {self.multiple_num_components:d}'

        return ret

#        TODO
#        star.cat_match = None
#        star.num_img_total = None
#        star.num_img_used = None
#        star.num_cat_pm = None
#        star.ra_mean_epoch = None
#        star.dec_mean_epoch = None
#        star.id_str = None
#        star.id_str_ucac2 = None


# --------------------------------------------------------------------------------
#    Bytes Format  Units   Label    Explanations
# --------------------------------------------------------------------------------
#    1-  4  I4     ---     HR       [1/9110]+ Harvard Revised Number
#                                     = Bright Star Number
#    5- 14  A10    ---     Name     Name, generally Bayer and/or Flamsteed name
#   15- 25  A11    ---     DM       Durchmusterung Identification (zone in
#                                     bytes 17-19)
#   26- 31  I6     ---     HD       [1/225300]? Henry Draper Catalog Number
#   32- 37  I6     ---     SAO      [1/258997]? SAO Catalog Number
#   38- 41  I4     ---     FK5      ? FK5 star Number
#       42  A1     ---     IRflag   [I] I if infrared source
#       43  A1     ---   r_IRflag  *[ ':] Coded reference for infrared source
#       44  A1     ---    Multiple *[AWDIRS] Double or multiple-star code
#   45- 49  A5     ---     ADS      Aitken's Double Star Catalog (ADS) designation
#   50- 51  A2     ---     ADScomp  ADS number components
#   52- 60  A9     ---     VarID    Variable star identification
#   61- 62  I2     h       RAh1900  ?Hours RA, equinox B1900, epoch 1900.0 (1)
#   63- 64  I2     min     RAm1900  ?Minutes RA, equinox B1900, epoch 1900.0 (1)
#   65- 68  F4.1   s       RAs1900  ?Seconds RA, equinox B1900, epoch 1900.0 (1)
#       69  A1     ---     DE-1900  ?Sign Dec, equinox B1900, epoch 1900.0 (1)
#   70- 71  I2     deg     DEd1900  ?Degrees Dec, equinox B1900, epoch 1900.0 (1)
#   72- 73  I2     arcmin  DEm1900  ?Minutes Dec, equinox B1900, epoch 1900.0 (1)
#   74- 75  I2     arcsec  DEs1900  ?Seconds Dec, equinox B1900, epoch 1900.0 (1)
#   76- 77  I2     h       RAh      ?Hours RA, equinox J2000, epoch 2000.0 (1)
#   78- 79  I2     min     RAm      ?Minutes RA, equinox J2000, epoch 2000.0 (1)
#   80- 83  F4.1   s       RAs      ?Seconds RA, equinox J2000, epoch 2000.0 (1)
#       84  A1     ---     DE-      ?Sign Dec, equinox J2000, epoch 2000.0 (1)
#   85- 86  I2     deg     DEd      ?Degrees Dec, equinox J2000, epoch 2000.0 (1)
#   87- 88  I2     arcmin  DEm      ?Minutes Dec, equinox J2000, epoch 2000.0 (1)
#   89- 90  I2     arcsec  DEs      ?Seconds Dec, equinox J2000, epoch 2000.0 (1)
#   91- 96  F6.2   deg     GLON     ?Galactic longitude (1)
#   97-102  F6.2   deg     GLAT     ?Galactic latitude (1)
#  103-107  F5.2   mag     Vmag     ?Visual magnitude (1)
#      108  A1     ---   n_Vmag    *[ HR] Visual magnitude code
#      109  A1     ---   u_Vmag     [ :?] Uncertainty flag on V
#  110-114  F5.2   mag     B-V      ? B-V color in the UBV system
#      115  A1     ---   u_B-V      [ :?] Uncertainty flag on B-V
#  116-120  F5.2   mag     U-B      ? U-B color in the UBV system
#      121  A1     ---   u_U-B      [ :?] Uncertainty flag on U-B
#  122-126  F5.2   mag     R-I      ? R-I   in system specified by n_R-I
#      127  A1     ---   n_R-I      [CE:?D] Code for R-I system (Cousin, Eggen)
#  128-147  A20    ---     SpType   Spectral type
#      148  A1     ---   n_SpType   [evt] Spectral type code
#  149-154  F6.3 arcsec/yr pmRA    *?Annual proper motion in RA J2000, FK5 system
#  155-160  F6.3 arcsec/yr pmDE     ?Annual proper motion in Dec J2000, FK5 system
#      161  A1     ---   n_Parallax [D] D indicates a dynamical parallax,
#                                     otherwise a trigonometric parallax
#  162-166  F5.3   arcsec  Parallax ? Trigonometric parallax (unless n_Parallax)
#  167-170  I4     km/s    RadVel   ? Heliocentric Radial Velocity
#  171-174  A4     ---   n_RadVel  *[V?SB123O ] Radial velocity comments
#  175-176  A2     ---   l_RotVel   [<=> ] Rotational velocity limit characters
#  177-179  I3     km/s    RotVel   ? Rotational velocity, v sin i
#      180  A1     ---   u_RotVel   [ :v] uncertainty and variability flag on
#                                     RotVel
#  181-184  F4.1   mag     Dmag     ? Magnitude difference of double,
#                                     or brightest multiple
#  185-190  F6.1   arcsec  Sep      ? Separation of components in Dmag
#                                     if occultation binary.
#  191-194  A4     ---     MultID   Identifications of components in Dmag
#  195-196  I2     ---     MultCnt  ? Number of components assigned to a multiple
#      197  A1     ---     NoteFlag [*] a star indicates that there is a note
#                                     (see file notes)
# --------------------------------------------------------------------------------
# Note (1): These fields are all blanks for stars removed from
#     the Bright Star Catalogue (see notes).
# Note on r_IRflag:
#   Blank if from NASA merged Infrared Catalogue, Schmitz et al., 1978;
#       ' if from Engles et al. 1982
#       : if uncertain identification
# Note on Multiple:
#       A = Astrometric binary
#       D = Duplicity discovered by occultation;
#       I = Innes, Southern Double Star Catalogue (1927)
#       R = Rossiter, Michigan Publ. 9, 1955
#       S = Duplicity discovered by speckle interferometry.
#       W = Worley (1978) update of the IDS;
# Note on n_Vmag:
#   blank = V on UBV Johnson system;
#       R = HR magnitudes reduced to the UBV system;
#       H = original HR magnitude.
# Note on pmRA:
#      As usually assumed, the proper motion in RA is the projected
#      motion (cos(DE).d(RA)/dt), i.e. the total proper motion is
#      sqrt(pmRA^2^+pmDE^2^)
# Note on n_RadVel:
#      V  = variable radial velocity;
#      V? = suspected variable radial velocity;
#      SB, SB1, SB2, SB3 = spectroscopic binaries,
#                          single, double or triple lined spectra;
#       O = orbital data available.
# --------------------------------------------------------------------------------

class YBSCStarCatalog(StarCatalog):
    def __init__(self,
                 dir: Optional[str | Path | FCPath] = None) -> None:
        """Create a YBSCStarCatalog.

        Parameters:
            dir: The path to the star catalog directory (may be a URL). Within
                 this directory should be the file ``catalog``.
        """

        super().__init__()

        if dir is None:
            self._dirname = FCPath(os.environ['YBSC_PATH'])
        else:
            self._dirname = FCPath(dir)

        self._stars = []

        with (self._dirname / 'catalog').open(mode='r') as fp:
            while True:
                record = fp.readline().rstrip()
                if record == '':
                    break
                record = record.ljust(197, ' ')
                if record[102:107].strip() == '':  # No VMAG
                    continue
                star = self._record_to_star(record)
                self._stars.append(star)

    def _find_stars(self,
                    ra_min: float,
                    ra_max: float,
                    dec_min: float,
                    dec_max: float,
                    vmag_min: Optional[float] = None,
                    vmag_max: Optional[float] = None,
                    full_result: bool = True,
                    **kwargs: Any) -> Iterator[YBSCStar]:

        # We do this here instead of as specific arguments because it works better
        # with mypy
        allow_double: bool = kwargs.pop('allow_double', False)

        for star in self._stars:
            if star.ra is None or star.dec is None:
                continue
            if not ra_min <= star.ra <= ra_max:
                continue
            if not dec_min <= star.dec <= dec_max:
                continue
            if star.vmag is not None:
                if vmag_min and star.vmag < vmag_min:
                    continue
                if vmag_max and star.vmag > vmag_max:
                    continue
            if not allow_double and star.multiple_star_code != ' ':
                continue

            if self.debug_level:
                print(star)
                print('-' * 80)

            yield star

    @staticmethod
    def _record_to_star(record: str) -> YBSCStar:

        star = YBSCStar()

        ###################
        # CATALOG NUMBERS #
        ###################

#    1-  4  I4     ---     HR       [1/9110]+ Harvard Revised Number
#                                     = Bright Star Number
#    5- 14  A10    ---     Name     Name, generally Bayer and/or Flamsteed name
#   15- 25  A11    ---     DM       Durchmusterung Identification (zone in
#                                     bytes 17-19)
#   26- 31  I6     ---     HD       [1/225300]? Henry Draper Catalog Number
#   32- 37  I6     ---     SAO      [1/258997]? SAO Catalog Number
#   38- 41  I4     ---     FK5      ? FK5 star Number

        star.unique_number = int(record[0:4].strip())
        star.name = record[4:14].strip()
        star.durchmusterung_id = record[14:25].strip()
        if record[25:31].strip() != '':
            star.draper_number = int(record[25:31].strip())
        if record[31:37].strip() != '':
            star.sao_number = int(record[31:37].strip())
        if record[37:41].strip() != '':
            star.fk5_number = int(record[37:41].strip())

        ################
        # SOURCE FLAGS #
        ################

#       42  A1     ---     IRflag   [I] I if infrared source
#       43  A1     ---   r_IRflag  *[ ':] Coded reference for infrared source
# Note on r_IRflag:
#   Blank if from NASA merged Infrared Catalogue, Schmitz et al., 1978;
#       ' if from Engles et al. 1982
#       : if uncertain identification
#       44  A1     ---    Multiple *[AWDIRS] Double or multiple-star code
# Note on Multiple:
#       A = Astrometric binary
#       D = Duplicity discovered by occultation;
#       I = Innes, Southern Double Star Catalogue (1927)
#       R = Rossiter, Michigan Publ. 9, 1955
#       S = Duplicity discovered by speckle interferometry.
#       W = Worley (1978) update of the IDS;
#   45- 49  A5     ---     ADS      Aitken's Double Star Catalog (ADS) designation
#   50- 51  A2     ---     ADScomp  ADS number components
#   52- 60  A9     ---     VarID    Variable star identification

        star.ir_source = (record[41] == 'I')
        if record[42] == ' ':
            star.ir_source_ref = YBSCStar.YBSC_IR_NASA
        elif record[42] == '\'':
            star.ir_source_ref = YBSCStar.YBSC_IR_ENGLES
        elif record[42] == ':':
            star.ir_source_ref = YBSCStar.YBSC_IR_UNCERTAIN

        star.multiple_star_code = record[43]
        if record[44:49].strip() != '':
            star.aitken_designation = record[44:49].strip()
        if record[49:51].strip() != '':
            star.ads_components = record[49:51].strip()
        if record[51:60].strip() != '':
            star.variable_star_id = record[51:60].strip()

        ###########
        # RA, DEC #
        ###########

#   76- 77  I2     h       RAh      ?Hours RA, equinox J2000, epoch 2000.0 (1)
#   78- 79  I2     min     RAm      ?Minutes RA, equinox J2000, epoch 2000.0 (1)
#   80- 83  F4.1   s       RAs      ?Seconds RA, equinox J2000, epoch 2000.0 (1)
#       84  A1     ---     DE-      ?Sign Dec, equinox J2000, epoch 2000.0 (1)
#   85- 86  I2     deg     DEd      ?Degrees Dec, equinox J2000, epoch 2000.0 (1)
#   87- 88  I2     arcmin  DEm      ?Minutes Dec, equinox J2000, epoch 2000.0 (1)
#   89- 90  I2     arcsec  DEs      ?Seconds Dec, equinox J2000, epoch 2000.0 (1)

        ra_hr = float(record[75:77])
        ra_min = float(record[77:79])
        ra_sec = float(record[79:83])
        dec_deg = float(record[83:86])
        dec_min = float(record[86:88])
        dec_sec = float(record[88:90])

        sign = 1
        if dec_deg < 0:
            dec_deg = -dec_deg
            sign = -1

        star.ra = np.radians((ra_hr/24. + ra_min/24./60 + ra_sec/24./60/60)*360)
        star.dec = np.radians(sign*(dec_deg + dec_min/60. + dec_sec/3600.))

        ########################
        # GALACTIC COORDINATES #
        ########################

#   91- 96  F6.2   deg     GLON     ?Galactic longitude (1)
#   97-102  F6.2   deg     GLAT     ?Galactic latitude (1)

        star.galactic_longitude = np.radians(float(record[90:96]))
        star.galactic_latitude = np.radians(float(record[96:102]))

        ##############
        # MAGNITUDES #
        ##############

#  103-107  F5.2   mag     Vmag     ?Visual magnitude (1)
#      108  A1     ---   n_Vmag    *[ HR] Visual magnitude code
# Note on n_Vmag:
#   blank = V on UBV Johnson system;
#       R = HR magnitudes reduced to the UBV system;
#       H = original HR magnitude.
#      109  A1     ---   u_Vmag     [ :?] Uncertainty flag on V
#  110-114  F5.2   mag     B-V      ? B-V color in the UBV system
#      115  A1     ---   u_B-V      [ :?] Uncertainty flag on B-V
#  116-120  F5.2   mag     U-B      ? U-B color in the UBV system
#      121  A1     ---   u_U-B      [ :?] Uncertainty flag on U-B
#  122-126  F5.2   mag     R-I      ? R-I   in system specified by n_R-I
#      127  A1     ---   n_R-I      [CE:?D] Code for R-I system (Cousin, Eggen)

        star.vmag = float(record[102:107])
        star.vmag_code = record[107]
        star.vmag_uncertainty_flag = record[108]
        if record[109:114].strip() != '':
            star.b_v = float(record[109:114])
        if record[115:120].strip() != '':
            star.u_b = float(record[115:120])
        if record[121:126].strip() != '':
            star.r_i = float(record[121:126])
        star.r_i_code = record[126]

        ##################
        # SPECTRAL CLASS #
        ##################

#  128-147  A20    ---     SpType   Spectral type
#      148  A1     ---   n_SpType   [evt] Spectral type code

        star.spectral_class = record[127:147].strip()
        star.spectral_class_code = record[147]

        #######################
        # MOTION AND PARALLAX #
        #######################

#  149-154  F6.3 arcsec/yr pmRA    *?Annual proper motion in RA J2000, FK5 system
# Note on pmRA:
#      As usually assumed, the proper motion in RA is the projected
#      motion (cos(DE).d(RA)/dt), i.e. the total proper motion is
#      sqrt(pmRA^2^+pmDE^2^)
#  155-160  F6.3 arcsec/yr pmDE     ?Annual proper motion in Dec J2000, FK5 system
#      161  A1     ---   n_Parallax [D] D indicates a dynamical parallax,
#                                     otherwise a trigonometric parallax
#  162-166  F5.3   arcsec  Parallax ? Trigonometric parallax (unless n_Parallax)
#  167-170  I4     km/s    RadVel   ? Heliocentric Radial Velocity
#  171-174  A4     ---   n_RadVel  *[V?SB123O ] Radial velocity comments
# Note on n_RadVel:
#      V  = variable radial velocity;
#      V? = suspected variable radial velocity;
#      SB, SB1, SB2, SB3 = spectroscopic binaries,
#                          single, double or triple lined spectra;
#       O = orbital data available.
#  175-176  A2     ---   l_RotVel   [<=> ] Rotational velocity limit characters
#  177-179  I3     km/s    RotVel   ? Rotational velocity, v sin i
#      180  A1     ---   u_RotVel   [ :v] uncertainty and variability flag on
#                                     RotVel

        star.pm_rac = float(record[148:154]) * AS_TO_RAD * YEAR_TO_SEC
        star.pm_ra = star.pm_rac / np.cos(cast(float, star.dec))
        star.pm_dec = float(record[154:160]) * AS_TO_RAD * YEAR_TO_SEC

        star.parallax_type = record[160]
        if record[161:166].strip() != '':
            star.parallax = float(record[161:166])
        if record[166:170].strip() != '':
            star.radial_velocity = float(record[166:170])
        star.radial_velocity_comments = record[170:174].strip()
        star.rotational_velocity_limit = record[174:176].strip()
        if record[176:179].strip() != '':
            star.rotational_velocity = float(record[176:179])
        star.rotational_velocity_uncertainty_flag = record[179:180]

#  181-184  F4.1   mag     Dmag     ? Magnitude difference of double,
#                                     or brightest multiple
#  185-190  F6.1   arcsec  Sep      ? Separation of components in Dmag
#                                     if occultation binary.
#  191-194  A4     ---     MultID   Identifications of components in Dmag
#  195-196  I2     ---     MultCnt  ? Number of components assigned to a multiple
#      197  A1     ---     NoteFlag [*] a star indicates that there is a note
#                                     (see file notes)

        if record[180:184].strip() != '':
            star.double_mag_diff = float(record[180:184])
        if record[184:190].strip() != '':
            star.double_mag_sep = float(record[184:190]) * AS_TO_RAD
        star.double_mag_components = record[190:194].strip()
        if record[194:196].strip() != '':
            star.multiple_num_components = int(record[194:196])

            ##################################################
            # COMPUTE SPECTRAL CLASS AND SURFACE TEMPERATURE #
            ##################################################

        sclass = star.spectral_class
        if sclass[0] == 'g':
            sclass = sclass[1:]
        sclass = sclass[0:2].strip()
        star.temperature = Star.temperature_from_sclass(sclass)

        return star
