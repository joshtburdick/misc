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
    private Map<String, Integer> columnNames = new HashMap<String, Integer>();

    /** The tableau. */
    private Tableau tableau;

    /**
     * Constructor.
     *
     * @param columnNames  names of the columns
     */
    public SimplexSolver(Vector<String> columnNames) {
        // number the columns
        int i = 0;
        for (String name : columnNames)
            this.columnNames.put(name, i++);


    }



}
