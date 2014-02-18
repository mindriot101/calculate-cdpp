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
