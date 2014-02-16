from astropy.io import fits as pyfits
from calculate_cdpp import cdpp
import os

def wd2jd(wd):
    jd_ref = 2453005.5
    return (wd / 86400.) + jd_ref

def test_calculation():
    base_path = os.path.dirname(__file__)
    data_path = os.path.join(base_path, 'data')

    data = pyfits.getdata(os.path.join(data_path, 'wasplc.fits'), 'photometry')
    hjd = wd2jd(data['tmid'].astype(float))
    flux = data['tamflux2'].astype(float)

    timescale_hours = 2
    median_cdpp, rms_cdpp = cdpp(hjd, flux, timescale_hours)

    assert 0.06 < median_cdpp < 0.08
    assert 0.08 < rms_cdpp < 0.10
