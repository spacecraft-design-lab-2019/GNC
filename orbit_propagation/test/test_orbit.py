import numpy as np
# from orbit_prop_py.orbit import *
import pytest

@pytest.mark.skip(reason="no way of currently testing this")
def test_get_orbit_pos():
    # test epoch
    epoch = '2013-12-14T14:18:37.00'
    # test TLE
    line1 = ('1 25635U 99008B   13348.59627062  .00000650  00000-0  16622-3 0  9860')
    line2 = ('2 25635 096.4421 173.2395 0141189 010.0389 029.8678 14.46831495780970')
    TLE = {'line1': line1, 'line2': line2}
    # Run function
    np.testing.assert_allclose(test_get_orbit_pos(TLE, epoch), (-5236.501106068836, 1137.472711241124, 4541.50244985641), atol=1e-6)

@pytest.mark.skip(reason="no way of currently testing this")
def test_get_orbit_magnetic():
    # test epoch
    epoch = '2013-12-14T14:18:37.00'
    # test TLE
    line1 = ('1 25635U 99008B   13348.59627062  .00000650  00000-0  16622-3 0  9860')
    line2 = ('2 25635 096.4421 173.2395 0141189 010.0389 029.8678 14.46831495780970')
    TLE = {'line1': line1, 'line2': line2}
    # Run function
    np.testing.assert_allclose(test_get_orbit_magnetic(TLE, epoch), 
        (-0.06298895551653533, 0.024098876110463806, 4.691724203241417, 
        5.228653116191708, -21.512815067758225,29.288854517493746, 26.25341004107369), atol=1e-6)