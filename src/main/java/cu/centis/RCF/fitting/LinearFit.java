/*
 * Copyright 2017 IAEA.
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

import java.util.Objects;
import org.apache.commons.math3.exception.NoDataException;

/**
 *
 * @author alex.vergara
 */
public class LinearFit {

    static class FLinear extends RobustFitter.MyParametricUnivariateFunction {

        @Override
        public double value(double t, double... parameters) {
            return parameters[0] + parameters[1] * t;
        }

        // Jacobian matrix of the above. In this case, this is just an array of
        // partial derivatives of the above function, with one element for each parameter.
        @Override
        public double[] gradient(double t, double... parameters) throws NoDataException {
            final double a = parameters[0];
            final double b = parameters[1];

            final double a1 = 1;
            final double b1 = t;

            return new double[]{a1, b1};
        }

        @Override
        public String getEquation() {
            return "y = a + b * x";
        }
    }

    public static class LinearFitter extends RobustFitter.MyAbstractCurveFitter {

        private LinearFitter(double[] xpoints, double[] ypoints, double[] weights) {
            this.xData = xpoints.clone();
            this.yData = ypoints.clone();
            if (Objects.isNull(weights)) {
                this.weights = new double[yData.length];
                for (int i = 0; i < yData.length; i++) {
                    this.weights[i] = 1;
                }
            } else {
                this.weights = weights.clone();
            }
            this.function = new FLinear();
        }

        public static LinearFitter create(double[] xpoints, double[] ypoints, double[] weights) {
            return new LinearFitter(xpoints, ypoints, weights);
        }

        @Override
        public void fit() {
            // Using default initialization
            double[] initialGuess = new double[]{0.0, 0.0};
            fit(initialGuess);
        }

        @Override
        public String getName() {
            return "Linear Fit";
        }

    }
}
