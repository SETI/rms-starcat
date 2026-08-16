=================
The Star Catalogs
=================

Three catalogs are supported. They differ enormously in size and in what they record, so
the right choice depends on the job.

.. list-table::
   :header-rows: 1
   :widths: 20 15 25 40

   * - Catalog
     - Stars
     - Data
     - Best for
   * - :class:`~starcat.SpiceStarCatalog`
     - varies
     - one ``.bdb`` kernel
     - pipelines that already use SPICE kernels
   * - :class:`~starcat.YBSCStarCatalog`
     - ~9,100
     - one 1.7 MB text file
     - naked-eye stars, with rich per-star detail
   * - :class:`~starcat.UCAC4StarCatalog`
     - ~113,000,000
     - 900 binary files, 8.6 GB
     - faint stars, astrometry, proper motion

Yale Bright Star Catalog
========================

The YBSC (Hoffleit & Warren 1991, 5th revised edition) covers essentially every star
visible to the naked eye, about 9,100 of them down to magnitude 6.5. What it lacks in
depth it makes up in breadth of detail: Bayer and Flamsteed names, Henry Draper, SAO and
FK5 numbers, UBV photometry, full spectral types, parallax, radial and rotational
velocity, and multiplicity.

The whole catalog is a single file named ``catalog``, obtainable from
http://tdc-www.harvard.edu/catalogs/bsc5.html. It is read into memory when the catalog
object is constructed, which takes a moment; each search after that is a scan over the
in-memory list. Records with no visual magnitude -- the entries that have been withdrawn
from the catalog -- are dropped while reading, which is why the object holds 9,096 stars
rather than the 9,110 lines in the file.

:class:`~starcat.YBSCStarCatalog` adds one option:

``allow_double`` (bool, default ``False``)
   Include stars that carry a double or multiple star code. **The default excludes about
   1,600 stars**, more than a sixth of the catalog, which is a surprise if you are
   counting. Pass ``True`` to see all of them.

.. code-block:: python

   from starcat import YBSCStarCatalog

   cat = YBSCStarCatalog()
   cat.count_stars()                     # 7519
   cat.count_stars(allow_double=True)    # 9096

UCAC4
=====

UCAC4 (Zacharias et al. 2013) holds about 113 million stars and is complete from the
brightest stars down to about magnitude R=16. Positional errors are 15 to 20 mas in the
10 to 14 magnitude range. Most stars have proper motions, and over 50 million carry
five-band photometry from APASS.

The data is distributed as 900 binary zone files in a ``u4b`` subdirectory, each covering
0.2 degrees of declination and sorted by increasing right ascension. Only the zones that
overlap the declination range of a search are opened, and within a zone the read stops as
soon as the right ascension passes ``ra_max``. **A narrow declination range is therefore
much cheaper than a narrow right ascension range**: a search one degree tall touches five
zone files no matter how wide it is.

:class:`~starcat.UCAC4StarCatalog` adds these options:

``require_clean`` (bool, default ``True``)
   Skip stars flagged as possible streaks, or with high proper motion that could not be
   matched or that disagrees with PPMXL.

``allow_double`` (bool, default ``True``)
   Include stars flagged as a component of a double or as a blended image.

``allow_galaxy`` (bool, default ``False``)
   Include objects matched to the LEDA galaxy catalog or to the 2MASS extended source
   catalog. UCAC4 is not exclusively a star catalog.

``require_pm`` (bool, default ``True``)
   Return only stars that have a measured proper motion.

``return_everything`` (bool, default ``False``)
   Override all four of the above to return every object in the box.

``optimize_ra`` (bool, default ``True``)
   Binary-search each zone file for the first record at or after ``ra_min`` instead of
   reading from the beginning. Turning it off is only useful when testing the search
   itself.

.. code-block:: python

   import numpy as np
   from starcat import UCAC4StarCatalog

   cat = UCAC4StarCatalog()

   # Clean stars with proper motion, the default
   n_good = cat.count_stars(dec_min=np.radians(0.), dec_max=np.radians(0.2))

   # Every object in the same zone, galaxies included
   n_all = cat.count_stars(dec_min=np.radians(0.), dec_max=np.radians(0.2),
                           return_everything=True)

A UCAC4 star carries two magnitudes: :attr:`~starcat.Star.vmag` is the
aperture magnitude, which the catalog considers the more robust of the two, and
:attr:`~starcat.UCAC4Star.vmag_model` is the fit model magnitude, which is more
accurate for well-behaved stars. The magnitude limits of a search apply to
:attr:`~starcat.Star.vmag`. Either may be ``None`` where the catalog recorded
no usable value.

UCAC4 records no spectral type. Where APASS supplies both a B and a V magnitude, this
package estimates :attr:`~starcat.Star.spectral_class` from the resulting B-V
color and derives :attr:`~starcat.Star.temperature` from that; where it does
not, both are ``None``. Treat these as estimates, not as catalog data.

SPICE catalogs
==============

NAIF distributes several star catalogs as SPICE type 1 kernels, among them ``hipparcos``,
``ppm``, and ``tycho2``. The file is named for the catalog with a ``.bdb`` extension, and
is loaded through CSPICE when the object is created:

.. code-block:: python

   from starcat import SpiceStarCatalog

   cat = SpiceStarCatalog('hipparcos')

A SPICE catalog records position, position uncertainty, magnitude, spectral class, and a
catalog number, and nothing else -- in particular **no proper motion**, so
:attr:`~starcat.Star.pm_ra` and :attr:`~starcat.Star.pm_dec` are
always ``None`` and :meth:`~starcat.Star.ra_dec_with_pm` returns the J2000
position unchanged. Some spectral classes in these catalogs end in an asterisk; the
conversions in :class:`~starcat.Star` strip it.

The kernel stays loaded for the life of the process. There are no catalog-specific search
options.

Choosing between them
=====================

- Working with naked-eye stars, or want spectral types, parallax, or multiplicity that
  someone else has curated? **YBSC.**
- Need faint stars, precise astrometry, or proper motion? **UCAC4**, if you can host the
  8.6 GB.
- Already have SPICE kernels in the pipeline and only need positions and magnitudes?
  **SPICE.**

Nothing stops you from using several at once; they share an interface precisely so that
the calling code does not have to care which one it was handed.
