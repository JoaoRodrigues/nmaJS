// Diagonalization Routines using LAPACK(E)
#include <stdio.h>
#include <stdlib.h>

#include <emscripten.h>
#include <lapacke.h>

// void ssyevr_(char *jobz, char *range, char *uplo, int *n, float *a, int *lda,
//              float *vl, float *vu, int *il, int *iu, float *abstol, int *m,
//              float *w, float *z, int *ldz, int *isuppz, float *work, int *lwork,
//              int *iwork, int *liwork, int *info);

void EMSCRIPTEN_KEEPALIVE diagonalize(float *coords, int natoms, int neig, float *evals, float *evecs);

void EMSCRIPTEN_KEEPALIVE diagonalize(float *mtx, int natoms, int neig, float *evals, float *evecs)
{

    // int i;  // debug

    // Parameters
    lapack_int info;

    lapack_int n = natoms;
    lapack_int lda = n;
    lapack_int ldz = neig;

    // printf("natoms = %d | neig = %d \n", natoms, neig);
    lapack_int il = 1;     // 1st eigval idx
    lapack_int iu = neig; // last idx.

    float abstol = -1.0; // default machine value.

    float vl = 0.0;
    float vu = 0.0;
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

    // printf("\n The total number of eigenvalues found:%2i\n", m);
    // for (i = 0; i < neig; i++)
    // {
    //     printf("eval %d = %f\n", i, evals[i]);
    // }
}
