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
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

/**
 * Test of SparseVector
 */
public class SparseVectorTest extends VectorTestAbstract {

    @Override
    protected void createPrimary() throws Exception {
        int n = Utilities.getInt(1, max);
        int m = Math.min(Utilities.getInt(max), n);
        x = new SparseVector(n);
        xd = Utilities.populate(x, m);
    }

    @Test
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

    @Test
    public void testBug27() {
        double[] tfVector = {0.0, 0.5, 0.0, 0.4, 0.0};
        DenseVector dense = new DenseVector(tfVector, false);
        SparseVector vectorTF = new SparseVector(dense);
        vectorTF.compact();

        assertTrue(vectorTF.getUsed() == 2); // vectorTF.getUsed() returns 5

        for (Iterator<VectorEntry> it = vectorTF.iterator(); it.hasNext();) {
            VectorEntry ve = it.next();
            int index = ve.index();
            double value = ve.get();
            assertTrue(tfVector[index] == value);
        }
    }

    /**
     * Unit test checking that the sparse vector does not end up ever using more
     * than "size" elements.
     */
    @Test
    public void testOverAllocation() {
        for (int d = 0; d < 10; d++) {
            SparseVector v = new SparseVector(d, 0);
            assertEquals(0, v.index.length);
            assertEquals(0, v.data.length);

            // Fill with non-zero elements.
            for (int i = 0; i < d; i++) {
                v.set(i, 1.0 + i);
            }

            assertEquals(d, v.index.length);
            assertEquals(d, v.data.length);
        }
    }

    @Test
    public void testGetRawIndex() {
        SparseVector vector = new SparseVector(Integer.MAX_VALUE);
        int[] index = vector.getRawIndex();
        assertTrue(index != null);
        assertTrue(index.length == 0);
        assertSame(index, vector.index);
        assertEquals(index.length, vector.getRawData().length);

        vector.set(2, 1.0);
        vector.set(1, 0.0);
        vector.set(4, 2.0);

        index = vector.getRawIndex();
        assertSame(index, vector.index);
        assertEquals(index.length, vector.getRawData().length);

        // In this case, the raw index is larger than the used, so the raw
        // indices have more entries than the other one.
        assertTrue(index.length > vector.getUsed());
        assertTrue(index.length > vector.getIndex().length);
    }

    @Test
    public void testGetRawData() {
        SparseVector vector = new SparseVector(Integer.MAX_VALUE);
        double[] data = vector.getRawData();
        assertTrue(data != null);
        assertTrue(data.length == 0);
        assertSame(data, vector.data);
        assertEquals(data.length, vector.getRawIndex().length);

        vector.set(2, 1.0);
        vector.set(1, 0.0);
        vector.set(4, 2.0);

        data = vector.getRawData();
        assertSame(data, vector.data);
        assertEquals(data.length, vector.getRawIndex().length);

        // In this case, the raw index is larger than the used, so the raw
        // indices have more entries than the other one.
        assertTrue(data.length > vector.getUsed());
        assertTrue(data.length > vector.getIndex().length);
    }
}
