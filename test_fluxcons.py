# -*- coding: utf-8 -*-

import pytest
import os
import glob

from drptestbones.backbone import consume_queue_directory
from drptestbones.diff import fits_osiris_allclose
from drptestbones.fetchdata import get_test_file

from drptestbones.checkFlux import compareFlux


def test_skyline(drf_queue):
    """Test FITS flux """
    consume_queue_directory(drf_queue)
    dark_sub_file = os.path.join(drf_queue, "dark_sub/s150508_c002080_Kbb_050.fits")
    reduced_file = os.path.join(drf_queue, "reduce/s150508_c002080_Kbb_050.fits")
    edges_dat = os.path.join(drf_queue, "edges.dat")

    darkFlux, redFlux = compareFlux(dark_sub_file,reduced_file,edges_dat,pixFromMax=6)
    fractFlux = darkFlux/redFlux

    fits_osiris_allclose(output_file, expected_file, edges_dat)
    assert ()

    
