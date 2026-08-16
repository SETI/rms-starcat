==================
Overview and Setup
==================

Purpose
=======

Different star catalogs answer different questions. The Yale Bright Star Catalog knows a
great deal about each of the few thousand stars visible to the naked eye; UCAC4 knows
position and magnitude for a hundred million stars far too faint to see; the SPICE
catalogs are what a NAIF-based pipeline already has on disk. Each is distributed in its
own format, with its own field layout, units, and conventions for missing data.

This package puts one interface in front of all three. You ask a catalog for the stars
within a box on the sky, optionally within a magnitude range, and you get back
:class:`~starcat.Star` objects whose common attributes -- position, magnitude,
proper motion, spectral class -- mean the same thing and use the same units no matter
which catalog produced them. Catalog-specific fields remain available on the subclasses.

.. mermaid::

   classDiagram
       StarCatalog <|-- SpiceStarCatalog
       StarCatalog <|-- YBSCStarCatalog
       StarCatalog <|-- UCAC4StarCatalog
       Star <|-- SpiceStar
       Star <|-- YBSCStar
       Star <|-- UCAC4Star
       StarCatalog ..> Star : yields

Installation
============

The package requires Python 3.10 or later:

.. code-block:: sh

   pip install rms-starcat

The catalog data is *not* included -- it is far too large -- and has to be obtained
separately. See :doc:`user_guide_catalogs` for what each catalog needs.

Locating the catalog data
=========================

Every catalog class takes an optional ``dir`` argument naming the directory that holds
its data. If you omit it, the directory is taken from an environment variable instead:

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Class
     - Environment variable
     - Expected contents
   * - :class:`~starcat.UCAC4StarCatalog`
     - ``UCAC4_PATH``
     - a ``u4b`` subdirectory of zone files
   * - :class:`~starcat.YBSCStarCatalog`
     - ``YBSC_PATH``
     - a file named ``catalog``
   * - :class:`~starcat.SpiceStarCatalog`
     - ``SPICE_PATH``
     - a ``Stars`` subdirectory of ``.bdb`` files

:class:`~starcat.SpiceStarCatalog` falls back to ``OOPS_RESOURCES``, in which case
it looks in the ``SPICE/Stars`` subdirectory, and raises ``RuntimeError`` if neither
variable is set. The other two raise ``KeyError`` if their variable is missing.

The path may be a local directory or a URL, because it is handled by
`rms-filecache <https://rms-filecache.readthedocs.io/en/latest/>`_. Remote data is
downloaded on demand and cached locally, so a catalog can be read straight from cloud
storage:

.. code-block:: python

   from filecache import FCPath
   from starcat import YBSCStarCatalog

   cat = YBSCStarCatalog(FCPath('gs://my-bucket/star-catalogs/YBSC'))

Units and conventions
=====================

- Right ascension, declination, and their uncertainties are in **radians**, at the J2000
  epoch. RA runs from 0 to 2*pi and DEC from -pi/2 to +pi/2.
- Proper motions are in **radians per second**, which is an awkward unit for a quantity
  normally quoted in milliarcseconds per year, but it is what you want when you multiply
  it by an elapsed time in seconds. The constants
  :data:`~starcat.starcatalog.MAS_TO_RAD` and :data:`~starcat.starcatalog.YEAR_TO_SEC`
  convert between the two.
- Magnitudes are visual magnitudes; smaller is brighter.
- Surface temperatures are in K.
- Any attribute may be ``None``, meaning the catalog had no value for that star. This is
  common: no SPICE catalog records proper motion, and most UCAC4 stars have no
  photometry from which a spectral class could be estimated.

A first search
==============

This finds Vega in the Yale Bright Star Catalog by searching a box a fifth of a degree
across:

.. code-block:: python

   import numpy as np
   from starcat import YBSCStarCatalog

   cat = YBSCStarCatalog()

   ra_vega = 279.2333
   dec_vega = 38.7836
   stars = list(cat.find_stars(ra_min=np.radians(ra_vega-0.1),
                               ra_max=np.radians(ra_vega+0.1),
                               dec_min=np.radians(dec_vega-0.1),
                               dec_max=np.radians(dec_vega+0.1)))

   print(len(stars))            # 1
   print(stars[0].name)         # 3Alp Lyr
   print(stars[0].vmag)         # 0.03

Printing a star gives a multi-line summary of everything the catalog knew about it,
which is the quickest way to see what a given catalog provides:

.. code-block:: python

   print(stars[0])

Where to go next
================

- :doc:`user_guide_searching` -- how the search box, the magnitude limits, and the
  performance options behave.
- :doc:`user_guide_catalogs` -- what each catalog contains, the options it adds, and
  where to obtain its data.
- :doc:`user_guide_stars` -- the attributes of a star, proper motion, and the spectral
  class conversions.
