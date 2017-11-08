/*
 * Copyright 2017 CENTIS.
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
package cu.centis.RCF.fitting;

import org.apache.commons.math3.exception.NoDataException;
import cu.centis.RCF.MathUtils;
import static cu.centis.RCF.fitting.RobustFitter.FWHM;
import static cu.centis.RCF.fitting.RobustFitter.FWTM;

/**
 *
 * @author alex.vergara
 */
public class GaussianFit {

    static class FGauss extends RobustFitter.MyParametricUnivariateFunction {

        @Override
        public double value(double t, double... parameters) {
            final double a = parameters[0];
            final double b = parameters[1];
            final double c = parameters[2];
            final double d = parameters[3];
            final double diff = t - c;
            final double i2d2 = 1 / (2 * d * d);
            return a + (b - a) * Math.exp(-diff * diff * i2d2);
        }

        // Jacobian matrix of the above. In this case, this is just an array of
        // partial derivatives of the above function, with one element for each parameter.
        @Override
        public double[] gradient(double t, double... parameters) throws NoDataException {
            final double a = parameters[0];
            final double b = parameters[1];
            final double c = parameters[2];
            final double d = parameters[3];

            final double diff = t - c;
            final double i2d2 = 1 / (2 * d * d);
            final double gauss = Math.exp(-diff * diff * i2d2);

            final double a1 = 1 - gauss;
            final double b1 = gauss;
            final double c1 = (b - a) * gauss * (diff / (d * d));
            final double d1 = (b - a) * gauss * (diff * diff / (d * d * d));

            return new double[]{a1, b1, c1, d1};
        }

        @Override
        public String getEquation() {
            return "y = a + (b-a)*exp(-(x-c)*(x-c)/(2*d*d))";
        }
    }

    public static class GaussianFitter extends RobustFitter.MyAbstractCurveFitter {

        private GaussianFitter(double[] xpoints, double[] ypoints) {
            this.xData = xpoints.clone();
            this.yData = ypoints.clone();
            this.weights = new double[yData.length];
            for (int i = 0; i < yData.length; i++) {
                this.weights[i] = 1;
            }
            this.function = new FGauss();
        }

        public static GaussianFitter create(double[] xpoints, double[] ypoints) {
            return new GaussianFitter(xpoints, ypoints);
        }

        @Override
        public synchronized void fit() {
            // Using default initialization
            double[] initialGuess = new double[]{0.0, MathUtils.Max(yData), MathUtils.Mean(xData), 1.0};
            fit(initialGuess);
        }

        @Override
        public String getName() {
            return "Gaussian Fit";
        }

        public synchronized double[] getResolution() {
            return new double[]{params[3] * FWHM, params[3] * FWTM}; 
        }

    }
}
