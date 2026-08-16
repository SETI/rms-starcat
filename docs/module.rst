``starcat`` Module
==================

.. automodule:: starcat
    :member-order: bysource
    :members:
    :undoc-members:
    :special-members:
    :show-inheritance:
    :exclude-members: __dict__, __hash__, __module__, __weakref__, __annotations__, __abstractmethods__

.. py:currentmodule:: starcat

Spectral class tables
---------------------

Both tables cover the main sequence classes from ``O5`` through ``M8``, including the
half class ``O9.5``. Values that the source did not give have been linearly
interpolated. The temperature difference between main sequence and giant stars is small
enough to ignore for these purposes.

The data is from Zombeck, M. V., *Handbook of Space Astronomy and Astrophysics*,
Cambridge University Press, 2nd ed., pp. 68-70, as transcribed at
http://www.vendian.org/mncharity/dir3/starcolor/details.html

.. py:data:: SCLASS_TO_B_MINUS_V
   :type: dict[str, float]

   Maps a spectral class, such as ``'G2'``, to its B-V color in magnitudes. Used by
   :meth:`Star.bmv_from_sclass` and, in reverse, by :meth:`Star.sclass_from_bv`.

.. py:data:: SCLASS_TO_SURFACE_TEMP
   :type: dict[str, float]

   Maps a spectral class, such as ``'G2'``, to its surface temperature in K. Used by
   :meth:`Star.temperature_from_sclass`. The M class values are from
   https://arxiv.org/abs/0903.3371 rather than from Zombeck.

Unit conversion constants
-------------------------

These are not exported from the top-level package; import them from
``starcat.starcatalog`` when converting the radians and radians-per-second used
throughout this package into more conventional units.

.. py:currentmodule:: starcat.starcatalog

.. py:data:: AS_TO_DEG
   :type: float

   Degrees per arcsecond.

.. py:data:: AS_TO_RAD
   :type: float

   Radians per arcsecond.

.. py:data:: MAS_TO_DEG
   :type: float

   Degrees per milliarcsecond.

.. py:data:: MAS_TO_RAD
   :type: float

   Radians per milliarcsecond. Divide a proper motion by this and by
   :data:`YEAR_TO_SEC` to get milliarcseconds per year.

.. py:data:: YEAR_TO_SEC
   :type: float

   Julian years per second, that is, the reciprocal of the number of seconds in a
   365.25-day year.

.. py:data:: TWOPI
   :type: float

   Two pi, the full range of right ascension in radians.

.. py:data:: HALFPI
   :type: float

   Pi over two, the declination of the north celestial pole in radians.
