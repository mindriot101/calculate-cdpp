import numpy as np
cimport numpy as np
import cython

cdef extern from "calculate_cdpp.h":
    void calculate_cdpp(double *hjd, double *flux, int N, double timescale_hours, double *median_cdpp, double *rms_cdpp)

@cython.boundscheck(False)
@cython.wraparound(False)
def cdpp(np.ndarray[double, ndim=1, mode="c"] hjd not None, 
        np.ndarray[double, ndim=1, mode="c"] flux not None, 
        double timescale_hours):
    '''
    Calculate the cdpp for a lightcurve.

    :param hjd:
        Time coordinate data in units of days

    :param flux:
        Flux measurements=

    :param timescale_hours:
        Hours to bin up to for the cdpp
    '''
    assert hjd.size == flux.size


    cdef double median_cdpp = 0.0
    cdef double rms_cdpp = 0.0

    cdef int size
    size = hjd.size

    calculate_cdpp(&hjd[0], &flux[0], size, timescale_hours, &median_cdpp, &rms_cdpp)

    return (median_cdpp, rms_cdpp)
