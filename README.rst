Polygonal Curve approximation
-----------------------------

Fast implementation of solutions to the min-num and min-e problems of digital
curve approximation.

Both approximate a planar curve given as a series of points by a subset of these
points with the minimum maximum euclidean distance between the original curve
and the linearly interpolated approximation. This distance is calculated
according to the tolerance zone criterion, described in `Alexander Kolesnikov's
thesis <http://cs.joensuu.fi/~koles/approximation/Ch3_1.html>`_

Both algorithms use the Ramer-Douglas-Peucker (RDP) algorithm as a subroutine.
This algorithm is fast and usually ends up with good approximations, but does
not guarantee optimality of the solution.

The min-e problem takes the length of the subset, i.e. the number of
approximating points and calculates the solution with the minimum maximum error.
Note that, due to the nature of the algorithm, the resulting subset length can
differ slightly from the required number of points.
The min-num problem takes an upper error bound and approximates the curve with
the minimum length subset meeting this error criterion.

The algorithms are implemented in C and accessible by a pythonic interface. The
min-num implementation offers similar functionality as the RDP implementations
by `sebleier <https://github.com/sebleier/RDP>`_ and `fhirschmann
<https://github.com/fhirschmann/rdp>`_, but can be orders of magnitudes faster
and more suitable for large curves. The present implementation uses a different
error measure than these implementations (tolerance zone vs infinite beam).


Installation
````````````
.. code:: bash

    pip install polyprox

Usage
`````

.. code:: python

    import numpy
    import polyprox

.. code:: python

    G = numpy.array([[0, 0], [0.9, 0], [1.1, 1.3], [2.5, 1.0], [2.2, 2.4]])
    polyprox.min_e(G, m=2)
    polyprox.min_num(G, epsilon=0.75)

Credits
```````

Kun Zhao: initial python implementation

License
```````

CSIRO Open Source Software License Agreement (variation of the BSD/MIT License)

References
``````````

| H. Imai and M. Iri. 1988. "Polygonal approximations of a curve - formulations and algorithms." In: Computational Morphology: A Computational Geometric Approach to the Analysis of Form: 71–86.

| W. S. Chan and F. Chin. 1996. "Approximation of polygonal curves with minimum number of line segments or minimum error." In: International Journal of Computational Geometry & Applications, 6 (01): 59–77

| Douglas, David H, and Thomas K Peucker. 1973. “Algorithms for the Reduction of the Number of Points Required to Represent a Digitized Line or Its Caricature.” In: Cartographica: The International Journal for Geographic Information and Geovisualization 10 (2): 112–122.
