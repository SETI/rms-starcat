==================
Working with Stars
==================

Every search yields :class:`~starcat.Star` objects, or rather objects of a
subclass: :class:`~starcat.SpiceStar`, :class:`~starcat.YBSCStar`, or
:class:`~starcat.UCAC4Star`. The attributes defined on the base class mean the same
thing and use the same units whichever catalog produced them.

The common attributes
=====================

- :attr:`~starcat.Star.unique_number` -- the catalog's own number for the
  star. It is unique within a catalog but says nothing across catalogs.
- :attr:`~starcat.Star.ra` and :attr:`~starcat.Star.dec` -- the
  J2000 position in radians, with uncertainties in
  :attr:`~starcat.Star.ra_sigma` and
  :attr:`~starcat.Star.dec_sigma`.
- :attr:`~starcat.Star.vmag` -- the visual magnitude, with its uncertainty in
  :attr:`~starcat.Star.vmag_sigma`.
- :attr:`~starcat.Star.pm_ra` and
  :attr:`~starcat.Star.pm_dec` -- proper motion in radians per second.
- :attr:`~starcat.Star.spectral_class` and
  :attr:`~starcat.Star.temperature` -- see below.

Attributes ending in ``rac`` rather than ``ra``, such as
:attr:`~starcat.Star.pm_rac`, are the right ascension quantity multiplied by
the cosine of the declination, which is how catalogs usually store them and what you want
for an angular rate on the sky rather than a rate of change of the coordinate.

**Any attribute may be** ``None``. Guard accordingly; a magnitude cut that reads
``star.vmag < 6`` will raise a ``TypeError`` on the first star whose magnitude the catalog
never measured.

Each subclass adds its own attributes on top -- see :doc:`user_guide_catalogs` and the
module reference. Accessing an attribute that does not exist raises ``AttributeError``,
as usual.

Proper motion
=============

Proper motion is stored in radians per second, which is not how it is usually quoted. To
read it in the conventional milliarcseconds per year, divide by the conversion constants:

.. code-block:: python

   from starcat.starcatalog import MAS_TO_RAD, YEAR_TO_SEC

   mas_per_year = star.pm_ra / MAS_TO_RAD / YEAR_TO_SEC

The reason for the choice is :meth:`~starcat.Star.ra_dec_with_pm`, which
moves a star to a different epoch given an elapsed time in seconds past J2000 -- the same
form as a SPICE ephemeris time:

.. code-block:: python

   import numpy as np

   # Where was this star at the start of 2025, 25 years after J2000?
   tdb = 25 * 365.25 * 86400.
   ra, dec = star.ra_dec_with_pm(tdb)
   print(np.degrees(ra), np.degrees(dec))

If the star has no proper motion, which is always the case for SPICE catalogs, the J2000
position is returned unchanged rather than an error; if it has no position at all, the
method returns ``(None, None)``. The correction is a simple linear extrapolation and
makes no attempt to handle the pole or to account for radial velocity.

Spectral class and temperature
==============================

What :attr:`~starcat.Star.spectral_class` contains depends on the catalog:

- **YBSC** gives the full spectral type as published, such as ``'A0Va'`` or
  ``'K0IIIbCN-0.5'``. The temperature is derived from the first two characters, after
  dropping a leading ``g`` for a giant.
- **SPICE** gives whatever the kernel holds, sometimes with a trailing ``*``.
- **UCAC4** publishes no spectral type at all. Where APASS supplies both a B and a V
  magnitude, this package estimates the class from the B-V color; otherwise both the
  class and the temperature are ``None``.

Three static methods on :class:`~starcat.Star` do the conversions, and are
useful on their own:

.. code-block:: python

   from starcat import Star

   Star.temperature_from_sclass('G2')     # 5780
   Star.bmv_from_sclass('G2')             # 0.63
   Star.sclass_from_bv(10.63, 10.)        # 'G2'

:meth:`~starcat.Star.temperature_from_sclass` and
:meth:`~starcat.Star.bmv_from_sclass` accept any capitalization, tolerate
surrounding whitespace and a trailing asterisk, and return ``None`` for a class they do
not recognize. :meth:`~starcat.Star.sclass_from_bv` returns the class whose
tabulated B-V is nearest, or ``None`` if the color falls outside the range of the table
(-0.32 to +2.00). Only the first two characters of a class are meaningful to these
methods, so pass ``'G2'`` and not ``'G2V'``.

The tables behind them, :data:`~starcat.SCLASS_TO_B_MINUS_V` and
:data:`~starcat.SCLASS_TO_SURFACE_TEMP`, are importable if you want to work with them
directly.

Converting to and from dictionaries
===================================

:meth:`~starcat.Star.to_dict` returns every attribute as a dictionary, and
:meth:`~starcat.Star.from_dict` sets attributes from one, which is convenient
for caching search results or moving them between processes:

.. code-block:: python

   import json

   records = [star.to_dict() for star in cat.find_stars(vmag_max=6.)]

   ...

   from starcat import YBSCStar

   stars = []
   for record in records:
       star = YBSCStar()
       star.from_dict(record)
       stars.append(star)

Two things to know. On a subclass the dictionary also contains the class constants, so a
:class:`~starcat.UCAC4Star` yields 66 entries rather than the 42 attributes it
actually carries; filter them out if that matters. And
:meth:`~starcat.Star.from_dict` sets whatever keys it is given, including
names that are not attributes of the class.

Printing a star
===============

``print(star)`` gives a multi-line summary of everything the catalog provided, with
angles shown in degrees and in sexagesimal notation and proper motions in milliarcseconds
per year. Missing values appear as ``None`` or ``N/A``. It is meant for reading, not for
parsing:

.. code-block:: text

   UNIQUE ID 7001 | RA 279.2345833° (18h36m56.300s) | DEC 38.7836111° (+038d47m1.000s)
   VMAG  0.030  | PM RA 259.135 mas/yr  | PM DEC 286.000 mas/yr
   TEMP 10800 | SCLASS A0Va
   Name "3Alp Lyr" | Durch "BD+38 3238" | Draper 172167 | SAO 67174 | FK5 699
   IR 1 Ref NASA | Multiple " " | Aitken 11510 None | Variable "Alp Lyr"
   SCLASS Code   | Galactic LON 67.44 LAT 19.24
   B-V 0.0 | U-B -0.01 | R-I -0.03
   Parallax TRIG 0.1230000 arcsec | RadVel -14.0 km/s V  | RotVel (v sin i) 15.0 km/s
   Double mag diff 10.40 Sep 62.80 arcsec Components AB # 5
