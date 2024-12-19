import unittest

import numpy as np

from filecache import FCPath

from starcat import SpiceStarCatalog


class Test_SpiceStarCatalog(unittest.TestCase):

    def runTest(self) -> None:
        cat = SpiceStarCatalog('hipparcos',
                               dir=FCPath('gs://rms-node-star-catalogs/SPICE/Stars',
                                          anonymous=True))

        num_all = cat.count_stars()
        self.assertEqual(num_all, 117955)

        num_vmag_lim = cat.count_stars(vmag_max=10)
        self.assertGreater(num_all, num_vmag_lim)

        # Compare slicing directions
        num_dec = 0
        for idec in range(20):
            num_dec += cat.count_stars(dec_min=np.radians(0.2*idec),
                                       dec_max=np.radians(0.2*(idec+1)),
                                       ra_min=np.radians(60), ra_max=np.radians(70))
        num_ra = 0
        for ira in range(10):
            num_ra += cat.count_stars(dec_min=np.radians(0.), dec_max=np.radians(4.),
                                      ra_min=np.radians(ira+60),
                                      ra_max=np.radians(ira+1+60))
        self.assertEqual(num_dec, num_ra)
