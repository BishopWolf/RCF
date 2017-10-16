/*
 * Copyright 2017 IAEA
 * Contract number TAL-NAHU20160831-001
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package cu.centis.RCF;

import ij.util.Tools;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Kurtosis;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;

/**
 * Enhanced version of StatUtils
 *
 * @author alex
 */
public class MathUtils {

    final static private double SQRT2 = Math.sqrt(2);

    /**
     * Private Constructor
     */
    private MathUtils() {
    }

    /**
     *
     * @param array the array values
     * @return the mean of the array values
     */
    public static double Mean(final double[] array)
            throws MathIllegalArgumentException {
        return StatUtils.mean(array);
    }

    /**
     *
     * @param array the array values
     * @return the mean of the array values
     */
    public static double Mean(final double[][] array)
            throws MathIllegalArgumentException {
        double tmean = Mean(array[0]);
        for (int i = 1; i < array.length; i++) {
            tmean += Mean(array[i]);
        }
        return tmean / array.length;
    }

    /**
     *
     * @param array the array values
     * @return the mean of the array values
     */
    public static double Median(final double[] array)
            throws MathIllegalArgumentException {
        return new Median().evaluate(array);
    }

    /**
     *
     * @param array the array values
     * @return the mode of the array values, if there are more than one mode
     * then they are returned in increasing order
     */
    public static double[] Mode(final double[] array)
            throws MathIllegalArgumentException {
        return StatUtils.mode(array);
    }

    /**
     *
     * @param array the array values
     * @return the min of the array values
     */
    public static double Min(final double[] array)
            throws MathIllegalArgumentException {
        return StatUtils.min(array);
    }

    /**
     * Modifies default min behavior so if second value is NAN return the first
     *
     * @param a the first value must not be NAN
     * @param b the second value to compare
     * @return if the second value is NAN return the first
     */
    public static double Min(double a, double b) {
        return Double.isNaN(b) || (b >= a) ? a : b;
    }

    /**
     * Modifies default min behavior so if second value is NAN return the first
     *
     * @param a the first value must not be NAN
     * @param b the second value to compare
     * @return if the second value is NAN return the first
     */
    public static float Min(float a, float b) {
        return Double.isNaN(b) || (b >= a) ? a : b;
    }

    /**
     * Modifies default min behavior so if second value is NAN return the first
     *
     * @param a the first value must not be NAN
     * @param b the second value to compare
     * @return if the second value is NAN return the first
     */
    public static int Min(int a, int b) {
        return Double.isNaN(b) || (b >= a) ? a : b;
    }

    /**
     *
     * @param array the array values
     * @return the min of the array values
     */
    public static double Min(final double[][] array)
            throws MathIllegalArgumentException {
        double tmin = Min(array[0]);
        for (int i = 1; i < array.length; i++) {
            tmin = Math.min(tmin, Min(array[i]));
        }
        return tmin;
    }

    /**
     *
     * @param array the array values
     * @return the min of the array values
     */
    public static float Min(final float[] array)
            throws MathIllegalArgumentException {
        return (float) Tools.getMinMax(array)[0];
    }

    /**
     *
     * @param array the array values
     * @return the min of the array values
     */
    public static float Min(final float[][] array)
            throws MathIllegalArgumentException {
        double[] a = Tools.getMinMax(array[0]);
        float tmin = (float) a[0];
        for (int i = 1; i < array.length; i++) {
            a = Tools.getMinMax(array[i]);
            tmin = (float) Math.min(tmin, a[0]);
        }
        return tmin;
    }

    /**
     *
     * @param array the array values
     * @return the min of the array absolute values
     */
    public static double ModularMin(final double[] array)
            throws MathIllegalArgumentException {
        double min = array[0];
        double mmin = Math.abs(min);
        for (double e : array) {
            if (Math.abs(e) < mmin && e != 0) {
                min = e;
                mmin = Math.abs(e);
            }
        }
        return min;
    }

    /**
     *
     * @param array the array values
     * @return the max of the array values
     */
    public static double Max(final double[] array)
            throws MathIllegalArgumentException {
        return StatUtils.max(array);
    }

    /**
     * Modifies default max behavior so if second value is NAN return the first
     *
     * @param a the first value must not be NAN
     * @param b the second value to compare
     * @return if the second value is NAN return the first
     */
    public static double Max(double a, double b) {
        return Double.isNaN(b) || (b <= a) ? a : b;
    }

    /**
     * Modifies default max behavior so if second value is NAN return the first
     *
     * @param a the first value must not be NAN
     * @param b the second value to compare
     * @return if the second value is NAN return the first
     */
    public static float Max(float a, float b) {
        return Double.isNaN(b) || (b <= a) ? a : b;
    }

    /**
     * Modifies default max behavior so if second value is NAN return the first
     *
     * @param a the first value must not be NAN
     * @param b the second value to compare
     * @return if the second value is NAN return the first
     */
    public static int Max(int a, int b) {
        return Double.isNaN(b) || (b <= a) ? a : b;
    }

    /**
     *
     * @param array the array values
     * @return the max of the array absolute values
     */
    public static double ModularMax(final double[] array)
            throws MathIllegalArgumentException {
        double max = array[0];
        double mmax = Math.abs(max);
        for (double e : array) {
            if (Math.abs(e) > mmax) {
                max = e;
                mmax = Math.abs(e);
            }
        }
        return max;
    }

    /**
     *
     * @param array the array values
     * @return the min of the array values
     */
    public static float Max(final float[] array)
            throws MathIllegalArgumentException {
        return (float) Tools.getMinMax(array)[1];
    }

    /**
     *
     * @param array the array values
     * @return the min of the array values
     */
    public static double Max(final double[][] array)
            throws MathIllegalArgumentException {
        double tmax = Max(array[0]);
        for (int i = 1; i < array.length; i++) {
            tmax = Math.min(tmax, Max(array[i]));
        }
        return tmax;
    }

    /**
     *
     * @param array the array values
     * @return the min of the array values
     */
    public static float Max(final float[][] array)
            throws MathIllegalArgumentException {
        double[] a = Tools.getMinMax(array[0]);
        float tmax = (float) a[1];
        for (int i = 1; i < array.length; i++) {
            a = Tools.getMinMax(array[i]);
            tmax = (float) Math.min(tmax, a[1]);
        }
        return tmax;
    }

    /**
     *
     * @param array the array values
     * @return the variance of the array values
     */
    public static double Variance(final double[] array)
            throws MathIllegalArgumentException {
        return StatUtils.variance(array);
    }

    /**
     *
     * @param array the array values
     * @return the skewness of the array values
     */
    public static double Skewness(final double[] array)
            throws MathIllegalArgumentException {
        return new Skewness().evaluate(array);
    }

    /**
     *
     * @param array the array values
     * @return the kurtosis of the array values
     */
    public static double Kurtosis(final double[] array)
            throws MathIllegalArgumentException {
        return new Kurtosis().evaluate(array);
    }

    /**
     *
     * @param array the array values
     * @return the standard deviation of the array values
     */
    public static double StdDev(final double[] array)
            throws MathIllegalArgumentException {
        return new StandardDeviation().evaluate(array);
    }

    /**
     *
     * @param array the array values
     * @return the sum of squares of the array values
     */
    public static double SumOfSquares(final double[] array)
            throws MathIllegalArgumentException {
        return StatUtils.sumSq(array);
    }

    /**
     *
     * @param array the array values
     * @return the mean of the array values
     */
    public static double GeometricMean(final double[] array)
            throws MathIllegalArgumentException {
        return StatUtils.geometricMean(array);
    }

    public static double[] Normalize(double[] array) {
        double norm2 = SumOfSquares(array);
        if (norm2 == 0) {
            return array;
        }
        double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i] / Math.sqrt(norm2);
        }
        return result;
    }

}
