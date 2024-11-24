package org.jburdick.exactsimplex;

import org.apache.commons.numbers.fraction.BigFraction;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

/**
 * A tableau for the simplex algorithm.
 */
public class Tableau {

    /** Names of columns. */
    private Map<String, Integer> columnNames;

    /**
     * The tableau, indexed by row, then column.
     *
     * We assume that the first row corresponds to the objective function.
     */
    public Map<Integer, Map<Integer, BigFraction>> t;

    /**
     * Constructor.
     *
     * @param columnNames  names of the columns
     */
    public Tableau(Vector<String> columnNames) {
    }





}
