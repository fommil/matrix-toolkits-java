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

import java.util.Iterator;
import java.util.List;

/**
 * An iterator over an array of iterable objects
 */
class SuperIterator<T extends Iterable<E>, E> implements
        Iterator<SuperIterator.SuperIteratorEntry> {

    private List<T> iterable;

    /**
     * Two iterators. We need the "next" iterator so that hasNext works properly
     * from one iterable to the next. Using a single iterator won't do
     */
    private Iterator<E> current, next;

    private int currentIndex = 0, nextIndex = 0;

    /**
     * Recyled entry returned from next()
     */
    private SuperIteratorEntry<E> entry;

    /**
     * Constructor for SuperIterator
     * 
     * @param iterable
     *            Iterable objects to iterate over
     */
    public SuperIterator(List<T> iterable) {
        this.iterable = iterable;
        entry = new SuperIteratorEntry<E>();

        // Try to be somewhat fault tolerant
        if (iterable.size() == 0) {
            current = new DummyIterator();
            next = new DummyIterator();
        } else {

            // This moves the next pointer to a non-empty iterable
            next = iterable.get(nextIndex).iterator();
            moveNext();

            // Then we move the current pointer in the same way
            current = iterable.get(currentIndex).iterator();
            moveCurrent();

            // Finally, move the next one step ahead if possible
            if (next.hasNext())
                next.next();
        }
    }

    private void moveNext() {
        while (nextIndex < iterable.size() - 1 && !next.hasNext())
            next = iterable.get(++nextIndex).iterator();
    }

    private void moveCurrent() {
        while (currentIndex < iterable.size() - 1 && !current.hasNext())
            current = iterable.get(++currentIndex).iterator();
    }

    public boolean hasNext() {
        return current.hasNext() || next.hasNext();
    }

    public SuperIteratorEntry<E> next() {
        // A wrapped object containing the relevant index and data
        entry.update(currentIndex, current.next());

        // Move current if necessary
        moveCurrent();

        // Move the next pointer
        moveNext();
        if (next.hasNext())
            next.next();

        return entry;
    }

    public void remove() {
        current.remove();
    }

    /**
     * Dummy iterator, for degenerate cases
     */
    private class DummyIterator implements Iterator<E> {

        public boolean hasNext() {
            return false;
        }

        public E next() {
            return null;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    /**
     * Entry returned from this superiterator
     */
    public static class SuperIteratorEntry<F> {

        /**
         * Index of the iterator which returned this
         */
        private int i;

        /**
         * Object returned
         */
        private F o;

        void update(int i, F o) {
            this.i = i;
            this.o = o;
        }

        public int index() {
            return i;
        }

        public F get() {
            return o;
        }

    }

}
