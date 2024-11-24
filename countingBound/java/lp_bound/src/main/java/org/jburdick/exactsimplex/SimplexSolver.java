package org.jburdick.exactsimplex;

import org.apache.commons.numbers.fraction.BigFraction;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

/**
 * Simplex algorithm.
 *
 * Following the description at
 * https://www.jeremykun.com/2014/12/01/linear-programming-and-the-simplex-algorithm/
 */
public class SimplexSolver {

    /** Names of columns. */
    private Map<String, Integer> columnNames;

    /**
     * The tableau, indexed by row, then column.
     */
    public Map<Integer, Map<Integer, BigFraction>> t;

    /**
     * Constructor.
     *
     * @param columnNames  names of the columns
     */
    public SimplexSolver(Vector<String> columnNames) {


    }




}
