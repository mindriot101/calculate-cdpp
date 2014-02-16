calculate-cdpp
==============

For the determination of the CDPP for a lightcurve.

Requirements
------------

* Cython
* Numpy
* cfitsio

Installation
------------

Possible options are

* download directory or clone git repository and run `python setup.py install`
* run `pip install git+git://github.com/mindriot101/calculate-cdpp.git`

Usage
-----

Given a lightcurve with hjd and flux measurements, and a desired timescale to bin over in hours, run

``` python
timescale_hours = 2.0
median_cdpp, rms_cdpp = calculate_cdpp.cdpp(hjd, flux, timescale_hours)
```

If you get errors such as invalid array types or contiguous memory then try

``` python
timescale_hours = 2.0
median_cdpp, rms_cdpp = calculate_cdpp.cdpp(hjd.astype(float), flux.astype(float), timescale_hours)
```

