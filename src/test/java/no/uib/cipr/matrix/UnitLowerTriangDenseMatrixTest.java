/*
 * Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 * 
 * This file is part of MTJ.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package no.uib.cipr.matrix;

/**
 * Test of UnitLowerTriangDenseMatrix
 */
public class UnitLowerTriangDenseMatrixTest extends UnitTriangMatrixTestAbstract {

    public UnitLowerTriangDenseMatrixTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        A = new UnitLowerTriangDenseMatrix(n);
        Ad = Utilities.unitLowerPopulate(A);
        Utilities.unitSet(Ad);
    }

  public void testMultUpper() {
    Matrix lu = new DenseMatrix(new double[][]{
        {-4.00,  6.00, 3.00},
        {1.00, -8.00,  5.00},
        {-0.50, -0.25,  0.75}
    });
    Matrix l = new UnitLowerTriangDenseMatrix(lu, false);
    Matrix u = new UpperTriangDenseMatrix(lu, false);

    Matrix e = new DenseMatrix(new double[][]{
        {-4,    6,    3},
        {-4,   -2,    8},
        {2,   -1,   -2}
    });

    Matrix out = l.mult(u, new DenseMatrix(3,3));
    MatrixTestAbstract.assertEquals(e, out);
  }

}
