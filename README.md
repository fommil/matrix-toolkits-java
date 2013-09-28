matrix-toolkits-java
====================

Java linear algebra library powered by BLAS and LAPACK

*MTJ* is designed to be used as a library for developing numerical applications, both for small and large scale computations. The library is based on [BLAS](http://www.netlib.org/blas) and [LAPACK](http://www.netlib.org/lapack) for its dense and structured sparse computations, and on the [Templates](http://www.netlib.org/templates) project for unstructured sparse operations.

MTJ uses the [netlib-java](https://github.com/fommil/netlib-java/) project as a backend,
which can be set up to use machine-optimised BLAS libraries for improved performance of
dense matrix operations.

History
=======

This project was originally written by Bjørn-Ove Heimsund, who has taken a step back due to other commitments. The group here is primarily concerned with keeping the library maintained, and fixing bugs as they are discovered. There is no road plan for future releases.

Installation
============

Releases are distributed on Maven central:

```xml
<dependency>
    <groupId>com.googlecode.matrix-toolkits-java</groupId>
    <artifactId>mtj</artifactId>
    <version>0.9.14</version>
</dependency>
```

Snapshots are distributed on Sonatype's Snapshot Repository:

```xml
<dependency>
  <groupId>com.googlecode.matrix-toolkits-java</groupId>
  <artifactId>mtj</artifactId>
  <version>1.0-SNAPSHOT</version>
</dependency>
```

Sparse Storage
==============

A variety of sparse matrix / vector storage classes are available:

* [`CompColMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/CompColMatrix.java)
* [`CompDiagMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/CompDiagMatrix.java)
* [`CompRowMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/CompRowMatrix.java)
* [`FlexCompColMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/FlexCompColMatrix.java)
* [`FlexCompRowMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/FlexCompRowMatrix.java)
* [`UnitLowerCompRowMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/UnitLowerCompRowMatrix.java)
* [`UpperCompRowMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/UpperCompRowMatrix.java)
* [`SparseVector`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/SparseVector.java)
* [`LinkedSparseMatrix`](blob/master/src/main/java/no/uib/cipr/matrix/sparse/LinkedSparseMatrix.java)

The `LinkedSparseMatrix` storage type is a novel storage type developed under this project. It maintains two tail links, one for the next matrix element by row order and another by column order. Lookups are kept into each row and column, making multiplication and transpose multiplication very fast.

The following charts compare the `LinkedSparseMatrix` against `DenseMatrix` for increasing matrix size (`n x n`) and number of non-zero elements, `m`. Lighter lines indicate larger `m`: varied from `10,000` to `100,000`. Solid lines are for dense matrix, dashed lines are the sparse matrix.

The following is time to initialise the matrix:

![init](http://i752.photobucket.com/albums/xx162/fommil/init_zpsec4b43b3.png)

The following is the memory consumption:

![mem](http://i752.photobucket.com/albums/xx162/fommil/mem_zpsc121c014.png)

The following is the time to perform a multiplication with a dense matrix and output into a dense matrix:

![mult](http://i752.photobucket.com/albums/xx162/fommil/mult_zps486a1af2.png)


Sparse Solvers
==============

MTJ provides [ARPACK](http://www.caam.rice.edu/software/ARPACK/) for very large symmetric matrices in [ArpackSym](blob/master/src/main/java/no/uib/cipr/matrix/sparse/ArpackSym.java) (see the example usage in [ArpackSymTest](blob/master/src/test/java/no/uib/cipr/matrix/sparse/ArpackSymTest.java)). ARPACK solves an arbitrary number of eigenvalues / eigenvectors.

In addition, several other solvers are available in the [`no.uib.cipr.matrix.sparse`](blob/master/src/test/java/no/uib/cipr/matrix/sparse) package. Users may wish to look at [Sparse Eigensolvers for Java](http://code.google.com/p/sparse-eigensolvers-java/) or consider implementing your own and letting us know about it (e.g. by using the [ARPACK](http://www.caam.rice.edu/software/ARPACK/) backend which comes with netlib-java).


Donations
=========

Please consider supporting the maintenance of this open source project with a donation:

[![Donate via Paypal](https://www.paypal.com/en_US/i/btn/btn_donateCC_LG.gif)](https://www.paypal.com/cgi-bin/webscr?cmd=_donations&business=B2HW5ATB8C3QW&lc=GB&item_name=mtj&currency_code=GBP&bn=PP%2dDonationsBF%3abtn_donateCC_LG%2egif%3aNonHosted)


Licence
=======

Copyright (C) 2003-2006 Bjørn-Ove Heimsund
Copyright (C) 2006-2013 Samuel Halliday

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see http://www.gnu.org/licenses/


Contributing
============

Contributors are encouraged to fork this repository and issue pull
requests. Contributors implicitly agree to assign an unrestricted licence
to Sam Halliday, but retain the copyright of their code (this means
we both have the freedom to update the licence for those contributions).
