import numpy as np
import unittest

from filecache import FCPath

from starcat import YBSCStarCatalog


class Test_YBSCStarCatalog(unittest.TestCase):

    def runTest(self) -> None:
        cat = YBSCStarCatalog(FCPath('gs://rms-node-star-catalogs/YBSC',
                                     anonymous=True))

        self.assertEqual(len(cat._stars), 9096)
        self.assertEqual(cat.count_stars(), 7519)
        self.assertEqual(cat.count_stars(allow_double=True), 9096)

        # Look up Vega
        ra_vega = 279.2333
        dec_vega = 38.7836

        vega_list = list(cat.find_stars(ra_min=np.radians(ra_vega-0.1),
                                        ra_max=np.radians(ra_vega+0.1),
                                        dec_min=np.radians(dec_vega-0.1),
                                        dec_max=np.radians(dec_vega+0.1)))
        self.assertEqual(len(vega_list), 1)

        vega = vega_list[0]
        self.assertEqual(vega.vmag, 0.03)
