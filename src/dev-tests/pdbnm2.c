// Calculates normal modes for a protein structure in PDB format
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <lapacke_utils.h>

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

float *calc_kirchhoff(float *coord, int n, float threshold, int *dd)
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
                mtx[si * d + sj] = -1.0; // u
                mtx[sj * d + si] = -1.0; // l
            }
        }
    }
    // for (i = 0; i < d; i++)
    // {
    //     for (j = 0; j < d; j++)
    //     {
    //         printf("mtx[%d, %d] (%d) = %4.1f\n", i, j, i * d + j, mtx[i * d + j]);
    //     }
    //     printf("\n");
    // }

    *dd = d;
    return mtx;
}

void diagonalize(float *mtx, int natoms, int neig, float *evals, float *evecs)
{

    int i;

    // Parameters
    lapack_int info;

    lapack_int n = natoms;
    lapack_int lda = n;
    lapack_int ldz = neig;

    printf("natoms = %d | neig = %d \n", natoms, neig);
    printf("The first 10 values of mtx are: ");
    for (i = 0; i < 10; i++)
    {
        printf("%f, ", mtx[i]);
    }
    printf("\n");
    lapack_int il = 1;    // 1st eigval idx
    lapack_int iu = neig; // last idx.

    float abstol = -1.0; // default machine value.

    float vl;
    float vu;
    lapack_int m;
    lapack_int isuppz[natoms];

    // Solve for eigenvalues/eigenvectors
    info = LAPACKE_ssyevr(
        LAPACK_ROW_MAJOR, 'V', 'I', 'U', n, mtx, lda,
        vl, vu, il, iu, abstol, &m, evals, evecs, ldz, isuppz);

    if (info > 0)
    {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }

    printf("\n The total number of eigenvalues found:%2i\n", m);
    for (i = 0; i < neig; i++)
    {
        printf("eval %d = %f\n", i, evals[i]);
    }
}

int main(int argc, char **argv)
{

    int n, d;
    float *mtx, *evals, *evecs;

    if (argc < 2)
    {
        printf("usage: pdbnm <pdbfile>\n");
        return 1;
    };

    float *coord = read_pdb(argv[1], &n);

    mtx = calc_kirchhoff(coord, n, 10.0, &d);

    int nmodes = strtol(argv[2], NULL, 10);

    evals = malloc(sizeof(float) * nmodes);
    evecs = malloc(sizeof(float) * nmodes * d);

    diagonalize(mtx, d, nmodes, evals, evecs);

    return 0;
}
