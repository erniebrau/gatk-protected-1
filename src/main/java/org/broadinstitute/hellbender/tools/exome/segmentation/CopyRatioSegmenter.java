package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.common.primitives.Doubles;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class CopyRatioSegmenter extends ScalarHMMSegmenter<Double> {
    private double logCoverageCauchyWidth;

    private static final double DEFAULT_INITIAL_CAUCHY_WIDTH = 0.1;
    private static final double MAX_REASONABLE_CAUCHY_WIDTH = 1.0;
    private static final double MIN_LOG_2_COPY_RATIO = -5.0;
    private static final double MAX_LOG_2_COPY_RATIO = 5.0;


    /**
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     */
    public CopyRatioSegmenter(final int initialNumStates, final ReadCountCollection rcc) {
        super(rcc.targets().stream().map(SimpleInterval::new).collect(Collectors.toList()), Doubles.asList(rcc.getColumn(0)),
                initialLog2CopyRatios(initialNumStates));
        Utils.validateArg(rcc.columnNames().size() == 1, "Only single-sample ReadCountCollection is supported.");
        logCoverageCauchyWidth = DEFAULT_INITIAL_CAUCHY_WIDTH;
    }

    public CopyRatioSegmenter(final int initialNumStates, final ReadCountCollection rcc, final double initialMemoryLength) {
        super(rcc.targets().stream().map(SimpleInterval::new).collect(Collectors.toList()), Doubles.asList(rcc.getColumn(0)),
                initialLog2CopyRatios(initialNumStates), initialMemoryLength);
        Utils.validateArg(rcc.columnNames().size() == 1, "Only single-sample ReadCountCollection is supported.");
        logCoverageCauchyWidth = DEFAULT_INITIAL_CAUCHY_WIDTH;
    }

    /**
     * evenly-spaced log-2 copy ratios
     * @param K the initial number of hidden states
     */
    private static List<Double> initialLog2CopyRatios(final int K) {
        final List<Double> result = new ArrayList<>();
        result.add(MIN_LOG_2_COPY_RATIO);       // 0
        result.add(ParamUtils.log2(0.5));       // 1
        result.add(0.0);   // 2
        IntStream.range(3, K).forEach(n -> result.add(ParamUtils.log2(n / 2.0)));
        return result;
    }

    @Override
    protected void relearnHiddenStateValues(final ExpectationStep eStep) { }

    public List<ModeledSegment> getModeledSegments() {
        final List<Pair<SimpleInterval, Double>> segmentation = findSegments();
        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        return segmentation.stream().map(pair ->
                new ModeledSegment(pair.getLeft(), tc.targetCount(pair.getLeft()), pair.getRight())).collect(Collectors.toList());
    }

    @Override
    protected ClusteringGenomicHMM<Double, Double> makeModel() {
        return new CopyRatioHiddenMarkovModel(getStates(), getMemoryLength(), logCoverageCauchyWidth);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        //relearn the Cauchy width of the emission distribution
        final Function<Double, Double> emissionLogLikelihood = width -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < numPositions(); position++) {
                for (int state = 0; state < numStates(); state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 : eStepPosterior
                            * CopyRatioHiddenMarkovModel.logEmissionProbability(data.get(position), getState(state), width);
                }
            }
            return logLikelihood;
        };

        logCoverageCauchyWidth = OptimizationUtils.argmax(emissionLogLikelihood, 0, MAX_REASONABLE_CAUCHY_WIDTH, logCoverageCauchyWidth,
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        logger.info("New coverage standard deviation learned: " + logCoverageCauchyWidth);
    }

    @Override
    protected double minHiddenStateValue() { return MIN_LOG_2_COPY_RATIO; }

    @Override
    protected double maxHiddenStateValue() { return  MAX_LOG_2_COPY_RATIO; }

    public double getLogCoverageCauchyWidth() {
        return logCoverageCauchyWidth;
    }

}
