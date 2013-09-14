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

/**
 * Array utilities. Complements <code>java.util.Arrays</code>
 * @deprecated java.utils.Arrays and Google Guava provide this functionality nowadays.
 */
@Deprecated
class Arrays {

    private Arrays() {
        // No need to instantiate
    }

    /**
     * Searches for a key in a sorted array, and returns an index to an element
     * which is greater than or equal key.
     * 
     * @param index
     *            Sorted array of integers
     * @param key
     *            Search for something equal or greater
     * @param begin
     *            Start posisiton in the index
     * @param end
     *            One past the end position in the index
     * @return end if nothing greater or equal was found, else an index
     *         satisfying the search criteria
     */
    public static int binarySearchGreater(int[] index, int key, int begin,
            int end) {
        return binarySearchInterval(index, key, begin, end, true);
    }

    /**
     * Searches for a key in a sorted array, and returns an index to an element
     * which is greater than or equal key.
     * 
     * @param index
     *            Sorted array of integers
     * @param key
     *            Search for something equal or greater
     * @return index.length if nothing greater or equal was found, else an index
     *         satisfying the search criteria
     */
    public static int binarySearchGreater(int[] index, int key) {
        return binarySearchInterval(index, key, 0, index.length, true);
    }

    /**
     * Searches for a key in a sorted array, and returns an index to an element
     * which is smaller than or equal key.
     * 
     * @param index
     *            Sorted array of integers
     * @param key
     *            Search for something equal or greater
     * @param begin
     *            Start posisiton in the index
     * @param end
     *            One past the end position in the index
     * @return begin-1 if nothing smaller or equal was found, else an index
     *         satisfying the search criteria
     */
    public static int binarySearchSmaller(int[] index, int key, int begin,
            int end) {
        return binarySearchInterval(index, key, begin, end, false);
    }

    /**
     * Searches for a key in a sorted array, and returns an index to an element
     * which is smaller than or equal key.
     * 
     * @param index
     *            Sorted array of integers
     * @param key
     *            Search for something equal or greater
     * @return -1 if nothing smaller or equal was found, else an index
     *         satisfying the search criteria
     */
    public static int binarySearchSmaller(int[] index, int key) {
        return binarySearchInterval(index, key, 0, index.length, false);
    }

    /**
     * Searches for a key in a subset of a sorted array.
     * 
     * @param index
     *            Sorted array of integers
     * @param key
     *            Key to search for
     * @param begin
     *            Start posisiton in the index
     * @param end
     *            One past the end position in the index
     * @return Integer index to key. -1 if not found
     */
    public static int binarySearch(int[] index, int key, int begin, int end) {
        return java.util.Arrays.binarySearch(index, begin, end, key);
    }

    private static int binarySearchInterval(int[] index, int key, int begin,
            int end, boolean greater) {

        // Zero length array?
        if (begin == end)
            if (greater)
                return end;
            else
                return begin - 1;

        end--; // Last index
        int mid = (end + begin) >> 1;

        // The usual binary search
        while (begin <= end) {
            mid = (end + begin) >> 1;

            if (index[mid] < key)
                begin = mid + 1;
            else if (index[mid] > key)
                end = mid - 1;
            else
                return mid;
        }

        // No direct match, but an inf/sup was found
        if ((greater && index[mid] >= key) || (!greater && index[mid] <= key))
            return mid;
        // No inf/sup, return at the end of the array
        else if (greater)
            return mid + 1; // One past end
        else
            return mid - 1; // One before start
    }

    /**
     * Finds the number of repeated entries
     * 
     * @param num
     *            Maximum index value
     * @param ind
     *            Indices to check for repetitions
     * @return Array of length <code>num</code> with the number of repeated
     *         indices of <code>ind</code>
     */
    public static int[] bandwidth(int num, int[] ind) {
        int[] nz = new int[num];

        for (int i = 0; i < ind.length; ++i)
            nz[ind[i]]++;

        return nz;
    }

}
