// Calculates normal modes for a protein structure in PDB format
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Declarations
void ssyevr_(char *jobz, char *range, char *uplo, int *n, float *a, int *lda,
             float *vl, float *vu, int *il, int *iu, float *abstol, int *m,
             float *w, float *z, int *ldz, int *isuppz, float *work, int *lwork,
             int *iwork, int *liwork, int *info);

/* Auxiliary routine: printing a matrix */
void print_matrix(char *desc, int m, int n, float *a, int lda)
{
    int i, j;
    printf("\n %s\n", desc);
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            printf(" %6.2f", a[i + j * lda]);
        printf("\n");
    }
}

// Main

float *read_pdb(char *fname, int *n)
{

    int i, j, s;

    FILE *fp;
    char line[80];
    int ncoord = 0;

    // What we want to read
    char record[5]; // extra char for null terminator
    char xyz[9];

    // Array of Coordinates
    int maxcoord = 1000;
    float *coord;
    coord = malloc(maxcoord * sizeof *coord); // init alloc

    fp = fopen(fname, "r");
    while (fgets(line, 80, fp))
    {

        memcpy(record, line, 4);
        record[4] = '\0';

        if (strcmp(record, "ATOM") != 0)
        {
            continue;
        }

        if (ncoord >= maxcoord)
        {
            maxcoord += 500; // grow array
            coord = realloc(coord, maxcoord * sizeof *coord);
        };

        // extract coordinates
        for (i = 0; i < 3; i++)
        {
            s = 30 + 8 * i; // stride
            for (j = 0; j < 8; j++)
            {
                xyz[j] = line[s + j];
            }
            coord[ncoord] = strtof(xyz, NULL);
            memset(xyz, 0, 9);
            ncoord++;
        }
    }

    fclose(fp);

    *n = ncoord;

    return coord;
}

float* calc_kirchhoff(float *coord, int n, float threshold, int *dd)
{
    int i, j;

    // init matrix
    int d = n / 3;
    float *mtx;
    mtx = malloc(d * d * sizeof *mtx);
    printf("n = %d | mtx d = %d | dd = %d\n", n, d, d * d);

    int si, sj;
    float sqthreshold = threshold * threshold;
    float xi, yi, zi, xj, yj, zj, sqdij;

    for (i = 0; i < n - 3; i += 3)
    {

        xi = coord[i];
        yi = coord[i + 1];
        zi = coord[i + 2];

        for (j = i + 3; j < n - 2; j += 3)
        {

            xj = coord[j];
            yj = coord[j + 1];
            zj = coord[j + 2];

            sqdij = ((xi - xj) * (xi - xj) +
                     (yi - yj) * (yi - yj) +
                     (zi - zj) * (zi - zj));

            if (sqdij <= sqthreshold)
            {

                si = i / 3; // strides in flattened array
                sj = j / 3;

                mtx[si * d + si] += 1.0;
                mtx[sj * d + sj] += 1.0;
                mtx[si * d + sj] = -1.0;  // u
                mtx[sj * d + si] = -1.0; // l
            }
        }
    }
    // for (i=0; i < d; i++) {
    //     for (j=0; j < d; j++) {
    //         printf("mtx[%d, %d] (%d) = %4.1f\n", i, j, i*d+j, mtx[i*d + j]);
    //     }
    //     printf("\n");
    // }

    *dd = d;
    return mtx;
}

int diagmtx(float *a, int n, int nsele)
{

    int i;

    // Parameters
    char jobz = 'V';
    char uplo = 'U';
    char range = 'I'; // to allow il-iu
    int lda = n;
    int ldz = n;

    int il = 1;     // 1st eigv
    int iu = nsele; // last eigv

    float vl;
    float vu;
    int m;  // num eigen vals
    float *w;  // eigen vals
    float *z;  // eigen vecs
    int *isuppz;
    float *work;
    int lwork;
    int liwork;
    int *iwork;
    int info;

    float wkopt;
    int iwkopt;

    w = malloc(sizeof(float) * n);
    printf("ldz = %d | nsele = %d\n", ldz, nsele);
    z = malloc(sizeof(float) * ldz * nsele);
    isuppz = malloc(sizeof(float) * n);

    /* Negative abstol means using the default value */
    float abstol = -1.0;

    // get workspace information:
    lwork = -1;
    liwork = -1;
    ssyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork,
            &info);

    printf("info: %d\n", info); // should be 0

    // do the real work:
    lwork = (int)wkopt;
    liwork = iwkopt;
    work = malloc(sizeof(float) * lwork);
    iwork = malloc(sizeof(int) * lwork);

    ssyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork,
            &info);

    printf("info: %d\n", info); // should be 0

    // find zero eigenvalues
    int nzm = 0;
    int *nzi = malloc(sizeof(int) * m);
    for (i=0; i < m; i++) {
        if ((w[i] - 0.0f) > 1e-5) {
            nzm++;
            nzi[i] = 1;
        }
    }

    /* Print the number of eigenvalues found */
    printf("\n The total number of eigenvalues found: %2i\n", m);
    printf("\n The total number of non-zero eigenvalues found: %2i\n", nzm);
    /* Print eigenvalues */
    print_matrix("Selected eigenvalues", 1, m, w, 1);
    /* Print eigenvectors */
    //   print_matrix("Selected eigenvectors (stored columnwise)", n, m, z, ldz);

    /* Free workspace */
    free((void *)iwork);
    free((void *)work);
    free((void *)nzi);
    free((void *)isuppz);

    return 0;
}

int main(int argc, char **argv)
{

    int n;
    int dd;
    float *mtx;

    if (argc < 2)
    {
        printf("usage: pdbnm <pdbfile>\n");
        return 1;
    };

    float *coord = read_pdb(argv[1], &n);

    printf("n is %d\n", n);
    mtx = calc_kirchhoff(coord, n, 10.0, &dd);

    int nmodes = strtol(argv[2], NULL, 10);
    printf("dd is %d\n", dd);
    diagmtx(mtx, dd, nmodes);


    return 0;
}
