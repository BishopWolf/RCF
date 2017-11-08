/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cu.centis.RCF.fitting;

import cu.centis.RCF.MathUtils;
import static cu.centis.RCF.fitting.RobustFitter.FWHM;
import static cu.centis.RCF.fitting.RobustFitter.FWTM;
import ij.plugin.filter.MaximumFinder;
import java.util.Arrays;
import org.apache.commons.math3.exception.NoDataException;

/**
 *
 * @author alex.vergara
 */
public class MultiGaussianFit {
    static class FMultiGauss extends RobustFitter.MyParametricUnivariateFunction {

        public int npeaks;

        public void setnpeaks(int lpeaks) {
            npeaks = lpeaks;
        }

        @Override
        public double value(double t, double... parameters) {
            double result = parameters[0];
            for (int i = 0; i < npeaks; i++) {
                final double b = parameters[3 * i + 1];
                final double c = parameters[3 * i + 2];
                final double d = parameters[3 * i + 3];
                final double diff = t - c;
                final double i2d2 = 1 / (2 * d * d);
                result += b * Math.exp(-diff * diff * i2d2);
            }
            return result;
        }

        // Jacobian matrix of the above. In this case, this is just an array of
        // partial derivatives of the above function, with one element for each parameter.
        @Override
        public double[] gradient(double t, double... parameters) throws NoDataException {
            double[] result = new double[3 * npeaks + 1];
            result[0] = 1;
            for (int i = 0; i < npeaks; i++) {
                final double b = parameters[3 * i + 1];
                final double c = parameters[3 * i + 2];
                final double d = parameters[3 * i + 3];

                final double diff = t - c;
                final double i2d2 = 1 / (2 * d * d);
                final double gauss = Math.exp(-diff * diff * i2d2);

                result[3 * i + 1] = gauss;
                result[3 * i + 2] = b * gauss * (diff / (d * d));
                result[3 * i + 3] = b * gauss * (diff * diff / (d * d * d));
            }
            return result;
        }

        @Override
        public String getEquation() {
            return String.format("y = a + sum[%1d, bi*exp(-(x-ci)*(x-ci)/(2*di*di))]", npeaks);
        }
    }

    public static class MultiGaussianFitter extends RobustFitter.MyAbstractCurveFitter {

        private MultiGaussianFitter(int npeaks, double[] xpoints, double[] ypoints) {
            this.xData = xpoints.clone();
            this.yData = ypoints.clone();
            this.weights = new double[yData.length];
            for (int i = 0; i < yData.length; i++) {
                this.weights[i] = 1;
            }
            this.function = new FMultiGauss();
            ((FMultiGauss) this.function).setnpeaks(npeaks);
        }

        public static MultiGaussianFitter create(int npeaks, double[] xpoints, double[] ypoints) {
            return new MultiGaussianFitter(npeaks, xpoints, ypoints);
        }

        @Override
        public synchronized void fit() {
            // Using default initialization
            int npeaks = ((FMultiGauss) this.function).npeaks;
            double[] initialGuess = new double[3 * npeaks + 1];
            initialGuess[0] = 0.0;
            double max = MathUtils.Max(yData);
            int[] peakpos = MaximumFinder.findMaxima(yData, 0.1 * max, false);
            Arrays.sort(peakpos);
            //TODO: handle overlapped peaks
            for (int i = 0; i < npeaks; i++) {
                initialGuess[3 * i + 1] = max;
                initialGuess[3 * i + 2] = xData[peakpos[i < peakpos.length ? i : 0]];
                initialGuess[3 * i + 3] = 1.0;
            }
            fit(initialGuess);
        }

        @Override
        public String getName() {
            return "Multi Gaussian Fit";
        }

        public synchronized double[] getResolution(int i) {
            return new double[]{params[3 * i + 3] * FWHM, params[3 * i + 3] * FWTM};
        }

    }
}
