package org.jburdick.exactsimplex;

import org.apache.commons.numbers.fraction.BigFraction;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

/**
 * A tableau for the simplex algorithm.
 */
public class Tableau {

    /** The tableau, indexed by row, then column.
     *
     * The last of each are indexed with -1.
     */
    public HashMap<Integer, HashMap<Integer, BigFraction>> t;

    /** Constructor. */
    public Tableau() {
        this.t = new HashMap<Integer, HashMap<Integer, BigFraction>>();
    }

    /** Adds one row. */
    public void addRow(HashMap<Integer, BigFraction> x) {


    }

}
