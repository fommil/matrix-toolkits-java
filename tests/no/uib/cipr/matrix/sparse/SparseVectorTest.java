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

package no.uib.cipr.matrix.sparse;

import java.util.Arrays;
import java.util.Iterator;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Utilities;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.VectorTestAbstract;

/**
 * Test of SparseVector
 */
public class SparseVectorTest extends VectorTestAbstract {

    public SparseVectorTest(String arg0) {
        super(arg0);
    }

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Math.min(Utilities.getInt(max), n);
        x = new SparseVector(n);
        xd = Utilities.populate(x, m);
    }

	public void testSparseVectorIndices() {
		/*
		 * MTJ subtlety in getIndex() for SparseVector. before calling
		 * getIndex(), you must call compact()... implementations may choose to
		 * do nothing in this call, but the Intel extended LAPACK
		 * implementations (and MTJ's SparseVector) require it. An alternative
		 * to vector.getIndex() is VectorMethods.getIndex(Vector) which will
		 * wrap this for you. It can take an arbitrary Vector and if it can be
		 * cast to a SparseVector will compact it and use its getIndex() method
		 * instead. (just so you're aware of this). Sam.
		 */

		// check that "infinite dimensions" doesn't use infinite memory
		SparseVector vector = new SparseVector(Integer.MAX_VALUE);
		int[] index = vector.getIndex();
		assert index != null;
		assert index.length == 0;

		// check that creating with double[] with zeros works
		double[] entries = new double[5];
		entries[0] = 0.0;
		entries[1] = 0.0;
		entries[2] = 1.0;
		entries[3] = 0.0;
		entries[4] = 2.0;
		Vector dense = new DenseVector(entries, false);
		vector = new SparseVector(dense);

		// NOTE: must compact before calling getIndex()!!!
		// vector.compact();
		index = vector.getIndex();
		assert index != null;
		assert index.length == 5 : "expected length of 5, but got "
				+ index.length + ", with elements " + Arrays.toString(index);
	}

	public void testBug27() {
		double[] tfVector = {0.0,0.5,0.0,0.4,0.0};
		DenseVector dense= new DenseVector(tfVector, false);
		SparseVector vectorTF =  new SparseVector(dense);
		vectorTF.compact();

		assertTrue(vectorTF.getUsed() == 2);  // vectorTF.getUsed() returns 5

		for (Iterator<VectorEntry> it = vectorTF.iterator();it.hasNext();) {
            VectorEntry ve= it.next();
            int index = ve.index();
            double value = ve.get();
            assertTrue(tfVector[index]== value);
        }
	}
}
