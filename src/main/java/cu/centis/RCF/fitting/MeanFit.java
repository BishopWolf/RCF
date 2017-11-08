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

/**
 *
 * @author alex.vergara
 */
public class MeanFit {

    static class FMean extends RobustFitter.MyParametricUnivariateFunction {

        @Override
        public double value(double t, double... parameters) {
            return parameters[0];
        }

        // Jacobian matrix of the above. In this case, this is just an array of
        // partial derivatives of the above function, with one element for each parameter.
        @Override
        public double[] gradient(double t, double... parameters) throws NoDataException {
            final double a = parameters[0];

            final double a1 = 1;

            return new double[]{a1};
        }

        @Override
        public String getEquation() {
            return "y = a";
        }
    }

    public static class MeanFitter extends RobustFitter.MyAbstractCurveFitter {

        private MeanFitter(double[] xpoints, double[] ypoints) {
            this.xData = xpoints.clone();
            this.yData = ypoints.clone();
            this.weights = new double[yData.length];
            for (int i = 0; i < yData.length; i++) {
                this.weights[i] = 1;
            }
            this.function = new FMean();
        }

        public static MeanFitter create(double[] xpoints, double[] ypoints) {
            return new MeanFitter(xpoints, ypoints);
        }

        @Override
        public synchronized void fit() {
            // Using default initialization
            double[] initialGuess = new double[]{0.0};
            fit(initialGuess);
        }

        @Override
        public String getName() {
            return "Mean Fit";
        }

    }
}
