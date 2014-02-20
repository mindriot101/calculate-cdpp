#include <Python.h>
#include <stdio.h>
#include <fitsio.h>
#include <stdlib.h>
#include <math.h>

#define SQ(x) ((x) * (x))

void fits_check(int status) {
    if (status) {
        fits_report_error(stderr, status);
        exit(status);
    }
}

void wd2jd(double *hjd, int N) {
    double jd_ref = 2453005.5;
    int i;
    for (i=0; i<N; ++i) {
        hjd[i] = (hjd[i] / 86400.0) + jd_ref;
    }
}

int compare(const void *a, const void *b) {
    const double *aa = (double*)a;
    const double *bb = (double*)b;

    if (*aa < *bb) {
        return -1;
    } else if (*aa > *bb) {
        return 1;
    } else {
        return 0;
    }
}

double mean(double *data, int N) {
    double total = 0;
    int i;
    for (i=0; i<N; ++i) {
        total += data[i];
    }
    return total / (double)N;
}

double _std(double *data, int N, double mean) {
    double total = 0;
    int i;
    for (i=0; i<N; ++i) {
        total += SQ(data[i] - mean);
    }
    return sqrt(total / (double)N);
}

double std(double *data, int N) {
    double mean_value = mean(data, N);
    return _std(data, N, mean_value);
}

void normalise(double *data, int N) {
    double mean_value = mean(data, N);
    int i;
    for (i=0; i<N; ++i) {
        data[i] /= mean_value;
        data[i] -= 1.0;
    }
}

void calculate_cdpp(double *hjd, double *flux, int N, double timescale_hours,
        double *median_cdpp, double *rms_cdpp) {

    *median_cdpp = 0;
    *rms_cdpp = 0;

    normalise(flux, N);

    double timescale_days = timescale_hours / 24.;
    double *cdpp_history = malloc(N * sizeof(double));

    int i, j;
    for (i=0; i<N; ++i) {
        double time_before = hjd[i] - timescale_days / 2.0;
        double time_after = hjd[i] + timescale_days / 2.0;

        long npts_in_bin = 0;
        double *within_bin = malloc(N * sizeof(double));


        for (j=0; j<N; ++j) {
            if ((hjd[j] >= time_before) && (hjd[j] <= time_after)) {
                within_bin[npts_in_bin++] = flux[j];
            }
        }
        within_bin = realloc(within_bin, npts_in_bin * sizeof(double));
        cdpp_history[i] = std(within_bin, npts_in_bin) / (double)sqrt(npts_in_bin);
        free(within_bin);
    }

    qsort(cdpp_history, N, sizeof(double), compare);

    *median_cdpp = cdpp_history[N/2];
    *rms_cdpp = std(cdpp_history, N);

    free(cdpp_history);

}

int main(int argc, const char *argv[]) {
    fitsfile *fptr;
    int status = 0;

    fits_open_file(&fptr, "data/wasplc.fits", READONLY, &status);
    fits_check(status);

    fits_movnam_hdu(fptr, BINARY_TBL, "PHOTOMETRY", 0, &status);
    fits_check(status);

    int hjd_col = 0;
    int flux_col = 0;
    fits_get_colnum(fptr, CASEINSEN, "TMID", &hjd_col, &status);
    fits_get_colnum(fptr, CASEINSEN, "TAMFLUX2", &flux_col, &status);
    fits_check(status);

    long nrows = 0;
    fits_get_num_rows(fptr, &nrows, &status);
    fits_check(status);

    printf("Reading %ld rows\n", nrows);

    double *hjd = malloc(nrows * sizeof(double));
    double *flux = malloc(nrows * sizeof(double));
    fits_read_col(fptr, TDOUBLE, hjd_col, 1, 1, nrows, NULL, hjd, NULL, &status);
    fits_read_col(fptr, TDOUBLE, flux_col, 1, 1, nrows, NULL, flux, NULL, &status);
    wd2jd(hjd, nrows);
    normalise(flux, nrows);
    fits_check(status);

    double timescale_hours = 2.;
    double median_cdpp = 0;
    double rms_cdpp = 0;
    calculate_cdpp(hjd, flux, nrows, timescale_hours, &median_cdpp, &rms_cdpp);

    printf("CDPP: %lf\n", median_cdpp);
    printf("RMS CDPP: %lf\n", rms_cdpp);


    free(hjd);
    free(flux);

    fits_close_file(fptr, &status);
    fits_check(status);
    
    return 0;
}
