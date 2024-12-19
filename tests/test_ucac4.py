import unittest

from filecache import FCPath
import numpy as np

from starcat import UCAC4StarCatalog


class Test_UCAC4StarCatalog(unittest.TestCase):

    def runTest(self) -> None:
        cat = UCAC4StarCatalog(FCPath('gs://rms-node-star-catalogs/UCAC4',
                                      anonymous=True))

        # Zone 1
        num_pm = cat.count_stars(require_clean=False, allow_double=True,
                                 allow_galaxy=True, require_pm=True,
                                 dec_min=np.radians(-90), dec_max=np.radians(-89.8))
        num_all = cat.count_stars(require_clean=False, allow_double=True,
                                  allow_galaxy=True, require_pm=False,
                                  dec_min=np.radians(-90), dec_max=np.radians(-89.8))
        self.assertEqual(num_all, 206)
        self.assertEqual(num_all-num_pm, 5)

        # Zone 451
        num_pm = cat.count_stars(require_clean=False, allow_double=True,
                                 allow_galaxy=True, require_pm=True,
                                 dec_min=np.radians(0.), dec_max=np.radians(0.2))
        num_all = cat.count_stars(require_clean=False, allow_double=True,
                                  allow_galaxy=True, require_pm=False,
                                  dec_min=np.radians(0.), dec_max=np.radians(0.2))
        self.assertEqual(num_all, 133410)
        self.assertEqual(num_all-num_pm, 6509)  # zone_stats says 6394??

        # Zone 900
        num_pm = cat.count_stars(require_clean=False, allow_double=True,
                                 allow_galaxy=True, require_pm=True,
                                 dec_min=np.radians(89.8), dec_max=np.radians(90))
        num_all = cat.count_stars(require_clean=False, allow_double=True,
                                  allow_galaxy=True, require_pm=False,
                                  dec_min=np.radians(89.8), dec_max=np.radians(90))
        self.assertEqual(num_all, 171)
        self.assertEqual(num_all-num_pm, 10)  # zone_stats says 9??

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
                                      ra_max=np.radians((ira+1)+60))
        self.assertEqual(num_dec, num_ra)

        # Compare optimized RA search with non-optimized
        for dec_idx in range(5):
            dec_min = np.radians(dec_idx*10-90.)
            dec_max = np.radians(dec_idx*10-89.9)
            for ra_min_idx in range(3):
                ra_min = np.radians(ra_min_idx * 10)
                ra_max = np.radians(ra_min + 10)
                num_opt = cat.count_stars(dec_min=dec_min, dec_max=dec_max,
                                          ra_min=ra_min, ra_max=ra_max)
                num_no_opt = cat.count_stars(dec_min=dec_min, dec_max=dec_max,
                                             ra_min=ra_min, ra_max=ra_max,
                                             optimize_ra=False)
                self.assertEqual(num_opt, num_no_opt)
