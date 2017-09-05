/*
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research
 * Organisation (CSIRO) - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the CSIRO Open Source Software License Agreement license.
 *
 * You should have received a copy of the CSIRO Open Source Software License
 * Agreement license with this file. If not, please write to:
 * kai.geissdoerfer@mailbox.tu-berlin.de
 */

#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Compare function for qsort */
int cmpf(const void * a, const void * b)
{
    double a_d = *(double*)a;
    double b_d = *(double*)b;
    if(a_d>b_d)
        return 1;
    else if(a_d<b_d)
        return -1;
    return 0;
}

/*
 * Infinite beam distance measure
 *
 * Calculates the minimum euclidean distance between the infinite line defined
 * by two points and a third point.
 *
 * x: x-coordinates of a curve
 * y: y-coordinates of a curve
 * start: Index of one of the points defining the line
 * end: Index of the other point defining the line
 * k: Index of third point
 *
 * returns: Euclidean distance between point k and the line
*/
double dist1(double * x, double * y, int start, int end, int k)
{
    double dx = (x[end]-x[start]);
    double dy = (y[end]-y[start]);
    double euc_dist = sqrt(dx*dx+dy*dy);


    return fabs((y[end]-y[start])*x[k]-(x[end]-x[start])*y[k]+x[end]*y[start]-y[end]*x[start])/euc_dist;
}

/*
 * Tolerance zone distance measure
 *
 * Calculates the minimum euclidean distance between a line segment defined by
 * two points and a third point.
 *
 * x: x-coordinates of a curve
 * y: y-coordinates of a curve
 * start: Index of one of the points defining the segment
 * end: Index of the other point defining the segment
 * k: Index of third point
 *
 * returns: Euclidean distance between point k and the line segment
*/
double dist2(double * x, double * y, int start, int end, int k)
{
    double gx = x[end] - x[start];
    double gy = y[end] - y[start];
    double g_l = gx*gx+gy*gy;

    double side_x = (x[k]-x[start]);
    double side_y = (y[k]-y[start]);
    double side_l = side_x*side_x+side_y*side_y;
    /*
     * If the two anchor points are relatively close, we can as well return the
     * distance to either of the anchor points. This avoids division by zero.
    */
    if(side_l > 1000.0*g_l)
        return sqrt(g_l);

    double p_g_l = (side_x*gx + side_y*gy)/g_l;

    if(p_g_l>1.0)
    {
        p_g_l = 1.0;
    }
    else if(p_g_l<0.0)
    {
        p_g_l = 0.0;
    }

    double p_g_x = x[start] + p_g_l * gx;
    double p_g_y = y[start] + p_g_l * gy;

    return sqrt((p_g_x-x[k])*(p_g_x-x[k]) + (p_g_y-y[k])*(p_g_y-y[k]));

}

/*
 * Ramer-Douglas-Peucker algorithm
 *
 * This algorithm solves the min-num problem. It approximates a digital curve by
 * a minimum length subset meeting a maximum error criterion, defined by
 * epsilon. The presented implementation uses a precalculated matrix, containing
 * the maximum error between any pair of points in the set (E) and the
 * corresponding index for this maximum error (I). Refer to the corresponding
 * Wikipedia page with pseudocode to understand the C implementation.
 *
 * E(optional): flattened matrix of maximum errors between every pair of nodes
 * I(optional): flattened matrix of maximum error indices between every pair of nodes
 * target: Buffer for indices of solution. Max length is number of nodes
 * epsilon: Maximum permitted error
 *
 * returns: Pointer to the next address after the last element written in target
*/

int * rdp(double * E, int * I, int * target, double * x, double * y, int start, int end, double epsilon)
{

    double d_max =  -1.0;
    int i_max = -1;

    if(E==NULL)
    {
        int k;
        /* If there's minimum of one hop between start and end */
        for(k=start+1;k<end;k++)
        {
            double d = dist2(x,y, start, end, k);
            if(d>d_max)
            {
                d_max = d;
                i_max = k;
            }
        }
    }
    else if(end-start>1)
    {
        /* Take precalculated error from flattened matrix */
        int index = (int)((end-2)*(end-1)/2 + start);
        d_max = E[index];
        i_max = I[index];
    }
    /* If the corresponding error exceeds epsilon, we need to add points in between */
    if (d_max > epsilon)
    {
        target = rdp(E, I, target, x, y, start, i_max, epsilon);
        /* Don't use the end point of the last recursion result */
        --target;
        target = rdp(E, I, target, x, y, i_max, end, epsilon);
    }
    else
    {
        *(target++) = start;
        *(target++) = end;
    }
    return target;
}

/*
 * Builds the error matrix between all pairs in the set.
 *
 * As can be seen easily, this matrix is less than a triangular matrix. In order
 * to save memory, the values are saved in a flattened array.
*/
double * build_matrix(double * x, double * y, int N, double * E, int * I)
{

    int i,j,k;

    for(i=0;i<N-2;i++)
    {
        for(j=i+2;j<N;j++)
        {

            double d_max = 0.0;
            int i_max = 0;
            /*
             * Iterate all vertices between the two points, for which the error
             * is calculated
            */
            for(k=i+1;k<j;k++)
            {
                double d = dist2(x,y, i, j, k);
                /* Search and store the maximum */
                if(d>d_max)
                {
                    d_max = d;
                    i_max = k;
                }
            }
            /*
             * We're filling something like a lower triangular matrix. In memory
             * we're operating on a flattened array by tricky indexing. The i-th
             * row and j-th element is in index i + \sum_{m=1}^{j-2} j
             * With gaussian sum formula this becomes (j-2)*(j-1)/2 + i
            */
            int index = (int)((j-2)*(j-1)/2 + i);
            E[index] = d_max;
            I[index] = i_max;
        }
    }

    return E;

}

/*
 * Implementation of the min-num problem solution
 *
 * This function is just the wrapper glueing above 'rdp' implementation to
 * python
*/
static PyObject* min_num_approximation(PyObject *dummy, PyObject *args)
{
    PyObject *arg1=NULL, *arg2=NULL;
    PyObject *arr1=NULL, *arr2=NULL;

    double epsilon;

    /* Parse input arrays */
    if (!PyArg_ParseTuple(args, "OOd", &arg1, &arg2, &epsilon)) return NULL;

    arr1 = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr1 == NULL) return NULL;
    arr2 = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr2 == NULL) goto fail;


    npy_intp * nx1 = PyArray_DIMS(arr1);
    int N = nx1[0];

    /* Convert input arrays to C arrays */
    double * dptrx = (double *)PyArray_DATA(arr1);
    double * dptry = (double *)PyArray_DATA(arr2);


    /* This will store the solution indices */
    int * P = malloc(N * sizeof(int));
    int * last = rdp(NULL, NULL, P, dptrx, dptry, 0, N-1, epsilon);
    int n_samples = 0;
    n_samples = (int) (last-P);

    npy_intp out_dims[1] = {n_samples};
    PyObject *narray = PyArray_SimpleNewFromData(1, out_dims, NPY_INT, P);

    /* Tell numpy it has to free the data */
    PyArray_ENABLEFLAGS((PyArrayObject*)narray, NPY_ARRAY_OWNDATA);

    Py_DECREF(arr1);
    Py_DECREF(arr2);

    Py_INCREF(Py_None);
    return narray;

 fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    return NULL;
}

/*
 * Implementation of the min-e problem solution
 *
 * Refer to 'Approximation of polygonal curves with minimum number of line
 * segments or minimum error' by CHAN and CHIN IJCGA Vol 6 1996 p. 59-77
 *
 * Takes the digital polygon defined by x and y and a number m of points.
 * Returns the minimum maximum error approximation of that curve with m samples
 * between start and end point.
*/
static PyObject* min_e_approximation(PyObject *dummy, PyObject *args)
{
    PyObject *arg1=NULL, *arg2=NULL;
    PyObject *arr1=NULL, *arr2=NULL;

    int m;

    // Parse input arrays
    if (!PyArg_ParseTuple(args, "OOi", &arg1, &arg2, &m)) return NULL;

    arr1 = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr1 == NULL) return NULL;
    arr2 = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr2 == NULL) goto fail;

    // Check input array lengths
    npy_intp * nx1 = PyArray_DIMS(arr1);
    npy_intp * nx2 = PyArray_DIMS(arr1);
    if(nx1[0]!=nx2[0])
        goto fail;


    int N = nx1[0];

    // Convert input arrays to C arrays
    double * dptrx = (double *)PyArray_DATA(arr1);
    double * dptry = (double *)PyArray_DATA(arr2);

    double * E = malloc((N-2)*(N-1)/2 * sizeof(double));
    int * I = malloc((N-2)*(N-1)/2 * sizeof(int));


    build_matrix(dptrx, dptry, N, E, I);

    /*
    The key idea here is, that we have already calculated all possible maximum
    errors (the maximum error corresponds to max error between two points).
    We now try those errors to see, by how many points they can be approximated
    with the min-num algorithm.
    If the number for a certain error is larger than we allow, we know, that we
    need to allow for a larger error and vice versa. This motivates below
    implementation of a binary search.
    */

    /* Copy flattened error matrix to new array of same size */
    double * E_sorted = malloc((N-2)*(N-1)/2 * sizeof(double));
    memcpy(E_sorted, E, (N-2)*(N-1)/2 * sizeof(double));

    /* Sort in ascending order */
    qsort(E_sorted, (N-2)*(N-1)/2, sizeof(double), cmpf);

    /* This will store the solution indices */
    int * P = malloc(N * sizeof(int));

    /* Initialize variables for a binary search */
    int mid, low = 0, high=(N-2)*(N-1)/2-1;

    int n_samples = 0;

    while(low<=high)
    {
        mid = (low+high)/2;
        /* Check by how many points the error can be achieved */
        int * last = rdp(E, I, P, dptrx, dptry, 0, N-1, E_sorted[mid]);
        n_samples = (int) (last-P);

        /* Our target number? Great, we're done */
        if(n_samples == m +2)
        {
            break;
        }
        /* Smaller number? Lower the allowed maximum error */
        else if(n_samples<(m+2))
        {
            high = mid - 1;
        }
        /* Bigger number? Relax the error */
        else
        {
            low = mid+1;
        }
    }


    npy_intp out_dims[1] = {n_samples};
    PyObject *narray = PyArray_SimpleNewFromData(1, out_dims, NPY_INT, P);
    /* Tell numpy it has to free the data */
    PyArray_ENABLEFLAGS((PyArrayObject*)narray, NPY_ARRAY_OWNDATA);

    free(E);
    free(E_sorted);
    free(I);

    /* Dereference the input arrays */
    Py_DECREF(arr1);
    Py_DECREF(arr2);

    Py_INCREF(Py_None);
    return narray;

 fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    return NULL;
}

static struct PyMethodDef methods[] = {
    {"min_e_approximation", min_e_approximation, METH_VARARGS, "Approximate given digitalized curve with minimum error by given number of samples"},
    {"min_num_approximation", min_num_approximation, METH_VARARGS, "Approximate given digitalized curve with minimum error by given number of samples"},
    {NULL, NULL}
};

PyMODINIT_FUNC
initalgorithms (void)
{
    (void)Py_InitModule("algorithms", methods);
    import_array();
}
