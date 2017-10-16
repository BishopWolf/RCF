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
import cu.centis.RCF.MathUtils;

/**
 *
 * @author alex.vergara
 */
public class SineFit {

    static class FSine extends RobustFitter.MyParametricUnivariateFunction {

        @Override
        public double value(double t, double... parameters) {
            return parameters[0] + parameters[1] * Math.sin(parameters[2] * t + parameters[3]);
        }

        // Jacobian matrix of the above. In this case, this is just an array of
        // partial derivatives of the above function, with one element for each parameter.
        @Override
        public double[] gradient(double t, double... parameters) throws NoDataException {
            final double a = parameters[0];
            final double b = parameters[1];
            final double c = parameters[2];
            final double d = parameters[3];

            final double a1 = 1;
            final double b1 = Math.sin(c * t + d);
            final double c1 = b * t * Math.cos(c * t + d);
            final double d1 = b * Math.cos(c * t + d);

            return new double[]{a1, b1, c1, d1};
        }

        @Override
        public String getEquation() {
            return "y = a + b * sin(c * x + d)";
        }
    }

    public static class SineFitter extends RobustFitter.MyAbstractCurveFitter {

        private SineFitter(double[] xpoints, double[] ypoints, double[] weights) {
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
            this.function = new FSine();
        }

        public static SineFitter create(double[] xpoints, double[] ypoints, double[] weights) {
            return new SineFitter(xpoints, ypoints, weights);
        }

        @Override
        public void fit() {
            // Using default initialization
            double[] initialGuess = new double[]{
                yData[0],
                0.5 * (MathUtils.Max(yData) - MathUtils.Min(yData)),
                2 * Math.PI / yData.length,
                0.0
            };
            fit(initialGuess);
        }

        @Override
        public String getName() {
            return "Sine Fit";
        }

    }

}
