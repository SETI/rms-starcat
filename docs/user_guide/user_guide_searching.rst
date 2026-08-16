=================
Searching the Sky
=================

Every catalog is searched the same way, through
:meth:`~starcat.StarCatalog.find_stars` and
:meth:`~starcat.StarCatalog.count_stars`. Both take the same criteria; the
first yields stars and the second only counts them.

The search box
==============

A search covers the box bounded by ``ra_min``, ``ra_max``, ``dec_min``, and ``dec_max``,
all in radians. The defaults cover the whole sky. All four are clipped to the valid range,
so passing a declination of -100 degrees is the same as passing -90.

.. code-block:: python

   import numpy as np
   from starcat import UCAC4StarCatalog

   cat = UCAC4StarCatalog()

   # A one degree square centered on RA 60, DEC 20
   stars = list(cat.find_stars(ra_min=np.radians(59.5), ra_max=np.radians(60.5),
                               dec_min=np.radians(19.5), dec_max=np.radians(20.5)))

Whether a star lying exactly on a boundary is included differs among the catalogs, so do
not depend on it. UCAC4 treats the box as closed at the bottom and open at the top; YBSC
includes both edges. If it matters, widen the box by a hair.

Wrapping across RA 0 and the poles
==================================

A box that crosses RA 0 is expressed by giving a ``ra_min`` *larger* than ``ra_max``. The
search is split internally into the two boxes that do not wrap:

.. code-block:: python

   # Ten degrees of RA centered on 0, i.e. 355 to 360 and 0 to 5
   stars = list(cat.find_stars(ra_min=np.radians(355.), ra_max=np.radians(5.),
                               dec_min=np.radians(19.5), dec_max=np.radians(20.5)))

The same is true of declination: a ``dec_min`` larger than ``dec_max`` selects the two
polar caps rather than the band between them. Giving both in the wrong order searches
four boxes.

Because a wrapped search is really several searches run one after another, the stars
come back grouped by sub-box. They are *not* in RA or DEC order across the seam, and in
the four-box case they are not in any single order at all. Sort the results yourself if
order matters.

Magnitude limits
================

``vmag_min`` and ``vmag_max`` restrict the visual magnitude. Remember that smaller
magnitudes are brighter, so ``vmag_max`` selects the *faintest* star you are willing to
accept:

.. code-block:: python

   # Everything brighter than 10th magnitude
   stars = list(cat.find_stars(vmag_max=10., ...))

   # A magnitude band
   stars = list(cat.find_stars(vmag_min=8., vmag_max=10., ...))

A limit of ``None`` means no limit; a limit of ``0.`` is a real limit, and excludes stars
on the wrong side of magnitude zero. A UCAC4 star whose magnitude is unknown fails any
limit you give and is dropped.

Counting rather than listing
============================

:meth:`~starcat.StarCatalog.count_stars` takes exactly the same criteria and
returns an integer. Prefer it to ``len(list(find_stars(...)))``: it passes
``full_result=False`` internally, so the catalogs skip the work of filling in fields
nobody is going to read.

.. code-block:: python

   n = cat.count_stars(dec_min=np.radians(19.5), dec_max=np.radians(20.5),
                       ra_min=np.radians(59.5), ra_max=np.radians(60.5))

Partial results
===============

You can ask for the same shortcut yourself with ``full_result=False``. The fields needed
to *decide* whether a star matches are always filled in -- position, magnitude, proper
motion, and the catalog's own quality flags -- but the rest are left as ``None``. For
UCAC4 that means skipping the uncertainties, the image and epoch counts, the APASS
photometry, the catalog match flags, the identification strings, and the derived
spectral class and temperature.

.. code-block:: python

   # Enough to plot positions, and appreciably faster
   for star in cat.find_stars(full_result=False, vmag_max=12., ...):
       plot(star.ra, star.dec)

Results arrive lazily
=====================

:meth:`~starcat.StarCatalog.find_stars` is a generator. Nothing is read until
you iterate, and a loop that stops early stops the reading. Iterate directly rather than
building a list when the box is large:

.. code-block:: python

   for star in cat.find_stars(vmag_max=6.):
       if star.vmag is not None and star.vmag < 0.:
           print(star.id_str)

Watching a search
=================

Setting ``debug_level`` on the catalog object prints a line for each star as it is kept
or rejected, along with the reason. UCAC4 also has a level 2, which additionally reports
the stars rejected on position and magnitude -- verbose, but the fastest way to find out
why a star you expected is missing.

.. code-block:: python

   cat.debug_level = 1
   cat.count_stars(...)
