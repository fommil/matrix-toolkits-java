matrix-toolkits-java
====================

Java linear algebra library powered by BLAS and LAPACK

*MTJ* is designed to be used as a library for developing numerical applications, both for small and large scale computations. The library is based on [BLAS](http://www.netlib.org/blas) and [LAPACK](http://www.netlib.org/lapack) for its dense and structured sparse computations, and on the [Templates](http://www.netlib.org/templates) project for unstructured sparse operations.

MTJ uses the [netlib-java](https://github.com/fommil/netlib-java/) project as a backend, which can be set up to use machine-optimised BLAS libraries for improved performance of dense matrix operations, falling back to a pure Java implementation. This ensures perfect portability, while allowing for improved performance in a production environment.

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

Sparse Solvers
==============

MTJ supports sparse matrix storage but does not provide solvers for sparse matrices. Have a look at [Sparse Eigensolvers for Java](http://code.google.com/p/sparse-eigensolvers-java/) or consider implementing your own and letting us know about it (e.g. by using the [ARPACK](http://www.caam.rice.edu/software/ARPACK/) backend which comes with netlib-java).


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
