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

import ij.IJ;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.util.Tools;
import java.awt.Color;
import java.util.Collection;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.linear.DiagonalMatrix;
import cu.centis.RCF.MathUtils;

/**
 *
 * @author alex.vergara
 */
public class RobustFitter {

    public static abstract class MyParametricUnivariateFunction implements ParametricUnivariateFunction {

        public abstract String getEquation();
    }

    public static abstract class MyAbstractCurveFitter extends AbstractCurveFitter {

        protected MyParametricUnivariateFunction function;
        protected double[] params;
        protected double[] xData, yData, weights;

        public void fit(double[] initialGuess) {
            this.params = initialGuess.clone();
            WeightedObservedPoints points = new WeightedObservedPoints();
            for (int i = 0; i < xData.length; i++) {
                points.add(new WeightedObservedPoint(weights[i], xData[i], yData[i]));
            }
            this.params = fit(points.toList());
        }

        public abstract void fit();

        public void setWeights(double[] lweights) {
            this.weights = lweights.clone();
        }

        public double f(double x) {
            return function.value(x, params);
        }

        public String getFormula() {
            return function.getEquation();
        }

        public abstract String getName();

        public double[] getParams() {
            return params;
        }

        public int getNumParams() {
            return params.length;
        }

        public double[] getXPoints() {
            return xData;
        }

        public double[] getYPoints() {
            return yData;
        }

        public double[] getResiduals() {
            double[] residuals = new double[xData.length];
            for (int i = 0; i < xData.length; i++) {
                residuals[i] = yData[i] - f(xData[i]);
            }
            return residuals;
        }

        private double getSumResidualsSqr() {
            double[] residuals = getResiduals();
            return MathUtils.Variance(residuals);
        }

        public double getRSquared() {
            double sumMeanDiffSqr = MathUtils.Variance(yData);
            double rSquared = 0.0;
            if (sumMeanDiffSqr > 0.0) {
                rSquared = 1.0 - getSumResidualsSqr() / sumMeanDiffSqr;
            }
            return rSquared;
        }

        /**
         * {@inheritDoc}
         */
        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> observations) {
            // Prepare least-squares problem.
            final int len = observations.size();
            final double[] target = new double[len];
            final double[] lweights = new double[len];

            int count = 0;
            for (WeightedObservedPoint obs : observations) {
                target[count] = obs.getY();
                lweights[count] = obs.getWeight();
                ++count;
            }

            final AbstractCurveFitter.TheoreticalValuesFunction model
                    = new AbstractCurveFitter.TheoreticalValuesFunction(function,
                            observations);

            // Create an optimizer for fitting the curve to the observed points.
            return new LeastSquaresBuilder().
                    maxEvaluations(Integer.MAX_VALUE).
                    maxIterations(Integer.MAX_VALUE).
                    start(params).
                    target(target).
                    weight(new DiagonalMatrix(lweights)).
                    model(model.getModelFunction(), model.getModelFunctionJacobian()).
                    build();
        }

        /**
         *
         * @param Title The title of the plot
         * @param xLabel label of x axis
         * @param yLabel label of y axis
         * @param eightBitCalibrationPlot resample image to 8 bit
         */
        public void plot(String Title, String xLabel, String yLabel, boolean eightBitCalibrationPlot) {
            String title = "".equals(Title) ? getFormula() : Title;
            if (getParams().length < getNumParams()) {
                Plot plot = new Plot(title, xLabel, yLabel, xData, yData);
                plot.setColor(Color.RED);
                plot.addLabel(0.02, 0.1, getName());
                plot.show();
                return;
            }
            int npoints = Math.min(Math.max(xData.length, 100), 1000);
            double[] a = Tools.getMinMax(xData);
            double xmin = a[0], xmax = a[1];
            if (eightBitCalibrationPlot) {
                npoints = 256;
                xmin = 0;
                xmax = 255;
            }
            a = Tools.getMinMax(yData);
            double ymin = a[0], ymax = a[1]; //y range of data points
            double[] px = new double[npoints];
            double[] py = new double[npoints];
            double inc = (xmax - xmin) / (npoints - 1);
            double tmp = xmin;
            for (int i = 0; i < npoints; i++) {
                px[i] = tmp;
                tmp += inc;
            }
            for (int i = 0; i < npoints; i++) {
                py[i] = f(px[i]);
            }
            a = Tools.getMinMax(py);
            double dataRange = ymax - ymin;
            ymin = Math.max(ymin - dataRange, Math.min(ymin, a[0])); //expand y range for curve, but not too much
            ymax = Math.min(ymax + dataRange, Math.max(ymax, a[1]));
            Plot plot = new Plot(title, xLabel, yLabel, px, py);
            plot.setLimits(xmin, xmax, ymin, ymax);
            plot.setColor(Color.RED);
            plot.addPoints(xData, yData, PlotWindow.CIRCLE);
            plot.setColor(Color.BLUE);

            StringBuilder legend = new StringBuilder(100);
            legend.append(getName()).append('\n');
            legend.append(getFormula()).append('\n');
            double[] p = getParams();
            int n = getNumParams();
            char pChar = 'a';
            for (int i = 0; i < n; i++) {
                legend.append(pChar).append(" = ").append(IJ.d2s(p[i], 5, 9)).append('\n');
                pChar++;
            }
            legend.append("R^2 = ").append(IJ.d2s(getRSquared(), 4)).append('\n');
            plot.addLabel(0.02, 0.1, legend.toString());
            plot.setColor(Color.BLUE);
            plot.show();
        }

        public void initializeIRLS() {
            fit();
            double[] Rx = getResiduals();
            for (int i = 0; i < Rx.length; i++) {
                Rx[i] = Math.abs(Rx[i]);
            }
            double[] xweights = MathUtils.Normalize(Rx);
            for (int i = 0; i < Rx.length; i++) {
                xweights[i] = 1 - xweights[i];
            }
            this.weights = xweights.clone();
        }

        public void runIRLS(int iterations) {
            int iter = 0;
            while (iter < iterations) {
                fit(params);
                if (getRSquared() > 0.9) {
                    break;
                }
                double[] Rx = getResiduals();
                for (int i = 0; i < Rx.length; i++) {
                    Rx[i] = Math.abs(Rx[i]);
                }
                double[] xweights = MathUtils.Normalize(Rx);
                for (int i = 0; i < Rx.length; i++) {
                    xweights[i] = 1 - xweights[i];
                }
                setWeights(xweights);
                ++iter;
            }
        }
    }

}
