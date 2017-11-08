/*
 * Copyright 2017 CENTIS
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
public class PolinomialFit {

    static class FPolinomial extends RobustFitter.MyParametricUnivariateFunction {

        private int order;

        public void setOrder(int lorder) {
            this.order = lorder;
        }

        public int getOrder() {
            return order;
        }

        @Override
        public double value(double t, double... parameters) {
            double result = 0.0;
            double var = 1.0;
            for (int i = 0; i <= order; i++) {
                result += parameters[i] * var;
                var *= t;
            }
            return result;
        }

        // Jacobian matrix of the above. In this case, this is just an array of
        // partial derivatives of the above function, with one element for each parameter.
        @Override
        public double[] gradient(double t, double... parameters) throws NoDataException {
            double[] result = new double[order + 1];
            double var = 1.0;
            for (int i = 0; i <= order; i++) {
                result[i] = var;
                var *= t;
            }
            return result;
        }

        @Override
        public String getEquation() {
            StringBuilder result = new StringBuilder();
            String var = " * x";
            result.append("y = ");
            for (int i = 0; i <= order; i++) {
                result.append("a").append(i);
                switch (i) {
                    case 0:
                        continue;
                    case 1:
                        result.append(var);
                        break;
                    default:
                        result.append(var).append(" ^ ").append(i);
                        break;
                }
            }
            return result.toString();
        }
    }

    public static class PolinomialFitter extends RobustFitter.MyAbstractCurveFitter {

        private PolinomialFitter(int order, double[] xpoints, double[] ypoints) {
            this.xData = xpoints.clone();
            this.yData = ypoints.clone();
                this.weights = new double[yData.length];
                for (int i = 0; i < yData.length; i++) {
                    this.weights[i] = 1;
                }
            this.function = new FPolinomial();
            ((FPolinomial) (this.function)).setOrder(order);
        }

        public static PolinomialFitter create(int order, double[] xpoints, double[] ypoints) {
            return new PolinomialFitter(order, xpoints, ypoints);
        }

        @Override
        public synchronized void fit() {
            // Using default initialization
            int order = ((FPolinomial) function).getOrder();
            double[] initialGuess = new double[order + 1];
            for (int i = 0; i <= order; i++) {
                initialGuess[i] = 1.0;
            }
            fit(initialGuess);
        }

        @Override
        public String getName() {
            return "Polinomial Fit";
        }

    }
}
