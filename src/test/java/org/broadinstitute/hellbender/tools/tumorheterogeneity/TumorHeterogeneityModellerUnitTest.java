package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Iterables;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.util.Log;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class TumorHeterogeneityModellerUnitTest extends BaseTest {
    private static final int RANDOM_SEED = 13;
    private static final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static final double CREDIBLE_INTERVAL_ALPHA = 0.95;

//    private static final File ACNV_SEG_FILE = new File("/home/slee/working/ipython/purity-ploidy/clonal_test_data/seed-1_trunc-frac-1.0_segments-1000_length-20/purity-0.7_total_segments.acnv.seg");
//    private static final File ACNV_SEG_FILE = new File("/home/slee/working/ipython/purity-ploidy/2_clones_test_data/seed-5_trunc-frac-0.5_segments-1000_length-20/purity-1.0_total_segments.acnv.seg");
    private static final List<File> ACNV_SEG_FILES = Arrays.asList(
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/0-SM-74NEG-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/1-SM-74P2T-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/2-SM-74P35-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/3-SM-74P3J-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/4-SM-74P3M-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/5-SM-74P3K-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/6-SM-74P51-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/7-SM-74P56-sim-final-edit.seg"), 
            new File("/home/slee/working/ipython/purity-ploidy/purity-series/8-SM-74P4M-sim-final-edit.seg"));

    @Test
    public void testRunMCMC() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        final int nMaxClonal = 6;
        final int nMax = 6;

        final int numPopulationsClonal = 2;
        final int numPopulations = 4;
        final int numCells = 400;

        final int numSamplesClonal = 50;
        final int numBurnInClonal = 25;
        final int numSamples = 100;
        final int numBurnIn = 50;

        final double concentrationPriorAlpha = 1.;
        final double concentrationPriorBeta = 1E6;
        final double variantSegmentFractionPriorAlpha = 4.;
        final double variantSegmentFractionPriorBeta = 4.;

        final PloidyState normalPloidyState = new PloidyState(1, 1);
        final Function<PloidyState, Double> ploidyPDF = ps -> Math.log(Math.pow(0.75, (ps.m() == 0 ? 1 : 0) + (ps.n() == 0 ? 1 : 0))) - 1000. * Math.log(Math.abs(normalPloidyState.m() - ps.m()) + Math.abs(normalPloidyState.n() - ps.n()));
//        final Function<PloidyState, Double> ploidyPDF = ps -> 0.;

        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMapClonal = new LinkedHashMap<>();
        for (int n = 0; n <= nMaxClonal; n++) {
            for (int m = 0; m <= n; m++) {
                if (m != 1 && n != 1) {
                    unnormalizedLogProbabilityMassFunctionMapClonal.put(new PloidyState(m, n), ploidyPDF.apply(new PloidyState(m, n)));
                }
            }
        }
        final PloidyStatePrior variantPloidyStatePriorClonal = new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMapClonal);

        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new LinkedHashMap<>();
        for (int n = 0; n <= nMax; n++) {
            for (int m = 0; m <= n; m++) {
                if (m != 1 && n != 1) {
                    unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(m, n), ploidyPDF.apply(new PloidyState(m, n)));
                }
            }
        }
        final PloidyStatePrior variantPloidyStatePrior = new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);

        for (final File ACNV_SEG_FILE : ACNV_SEG_FILES) {
            rng.setSeed(RANDOM_SEED);

            final String filePath = ACNV_SEG_FILE.getAbsolutePath().replace(".seg", "");
            final List<ACNVModeledSegment> allSegments = SegmentUtils.readACNVModeledSegmentFile(ACNV_SEG_FILE);

            final File resultClonalFile = new File(filePath + "-result-clonal.tsv");
            resultClonalFile.createNewFile();
            try (final FileWriter writerClonal = new FileWriter(resultClonalFile)) {
                output(writerClonal, logger, "#num segments all: " + allSegments.size());
                output(writerClonal, logger, System.getProperty("line.separator"));
                final List<ACNVModeledSegment> segments = filterSegments(allSegments, writerClonal, logger);
                output(writerClonal, logger, "#num segments after filtering: " + segments.size());
                output(writerClonal, logger, System.getProperty("line.separator"));

                final File mafCrFile = new File(filePath + "-maf-cr.tsv");
                mafCrFile.createNewFile();
                final TumorHeterogeneityData data = new TumorHeterogeneityData(segments);
                try (final FileWriter mafCrWriter = new FileWriter(mafCrFile)) {
                    for (double f = 0.; f <= 0.5; f += 0.01) {
                        for (double c = 1E-10; c <= 5; c += 0.05) {
                            final double maf = f;
                            final double cr = c;
                            final double density = IntStream.range(0, data.numSegments()).mapToDouble(i -> Math.exp(data.logDensity(i, cr, maf))).sum();
                            mafCrWriter.write(maf + "\t" + cr + "\t" + density + System.getProperty("line.separator"));
                        }
                    }
                }
                
                //run MCMC
                final TumorHeterogeneityModeller clonalModeller = new TumorHeterogeneityModeller(
                        data, normalPloidyState, variantPloidyStatePriorClonal,
                        concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta,
                        numPopulationsClonal, numCells, rng);
                clonalModeller.fitMCMC(numSamplesClonal, numBurnInClonal);
                outputModeller(ctx, segments, variantPloidyStatePriorClonal, numPopulationsClonal, numCells, numSamplesClonal, numBurnInClonal, clonalModeller, writerClonal, logger);

                final File resultFile = new File(filePath + "-result.tsv");
                resultFile.createNewFile();
                try (final FileWriter writer = new FileWriter(resultFile)) {
                    final double clonalConcentration = Iterables.getLast(clonalModeller.getConcentrationSamples());
                    final TumorHeterogeneityState.PopulationFractions initialPopulationFractions = new TumorHeterogeneityState.PopulationFractions(Collections.nCopies(numPopulations, 1. / numPopulations));
                    final TumorHeterogeneityState.PopulationIndicators initialPopulationIndicators = Iterables.getLast(clonalModeller.getPopulationIndicatorsSamples());
                    IntStream.range(0, numCells).filter(i -> initialPopulationIndicators.get(i) == 1).forEach(i -> initialPopulationIndicators.set(i, numPopulations - 1));
                    final TumorHeterogeneityState.VariantProfileCollection clonalVariantProfileCollection = Iterables.getLast(clonalModeller.getVariantProfileCollectionSamples());
                    final List<TumorHeterogeneityState.VariantProfile> initialVariantProfiles = new ArrayList<>();
                    initialVariantProfiles.addAll(clonalVariantProfileCollection);
                    initialVariantProfiles.add(1, TumorHeterogeneityModeller.initializeProfile(segments.size()));
                    initialVariantProfiles.add(1, TumorHeterogeneityModeller.initializeProfile(segments.size()));
                    final TumorHeterogeneityState.VariantProfileCollection initialVariantProfileCollection =
                            new TumorHeterogeneityState.VariantProfileCollection(initialVariantProfiles);
                    final TumorHeterogeneityModeller modeller = new TumorHeterogeneityModeller(
                            clonalConcentration, initialPopulationFractions, initialPopulationIndicators, initialVariantProfileCollection,
                            data, normalPloidyState, variantPloidyStatePrior,
                            concentrationPriorAlpha, concentrationPriorBeta, variantSegmentFractionPriorAlpha, variantSegmentFractionPriorBeta,
                            numPopulations, numCells, rng);
                    modeller.fitMCMC(numSamples, numBurnIn);
                    outputModeller(ctx, segments, variantPloidyStatePrior, numPopulations, numCells, numSamples, numBurnIn, modeller, writer, logger);
                }
            }
        }
    }
    
    private void output(final FileWriter writer, final Logger logger, final String output) {
        try {
            writer.write(output);
            logger.info(output);
        } catch(final IOException e) {
            throw new GATKException("Cannot output.");
        }
    }

    private List<ACNVModeledSegment> filterSegments(final List<ACNVModeledSegment> allSegments, final FileWriter writer, final Logger logger) {
        final double lengthPercentile = 0.1;
        final double credibleIntervalPercentile = 95.;

        final Percentile percentile = new Percentile();

        final double[] lengths = allSegments.stream().mapToDouble(s -> (double) s.getInterval().size()).toArray();
        final int lengthThreshold = (int) percentile.evaluate(lengths, lengthPercentile);
        output(writer, logger, "#length threshold: " + lengthThreshold);
        output(writer, logger, System.getProperty("line.separator"));

//        final double[] log2crCredibleIntervalSizes = allSegments.stream()
//                .mapToDouble(s -> s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_10))
//                .toArray();
//        final double log2crCredibleIntervalThreshold = percentile.evaluate(log2crCredibleIntervalSizes, credibleIntervalPercentile);
//        System.out.println("Log2CR credible-interval threshold: " + log2crCredibleIntervalThreshold);

        final double[] mafCredibleIntervalSizes = allSegments.stream()
                .filter(s -> !Double.isNaN(s.getMinorAlleleFractionPosteriorSummary().getCenter()))
                .mapToDouble(s -> s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_10))
                .toArray();
        final double mafCredibleIntervalThreshold = percentile.evaluate(mafCredibleIntervalSizes, credibleIntervalPercentile);
        output(writer, logger, "#MAF credible-interval threshold: " + mafCredibleIntervalThreshold);
        output(writer, logger, System.getProperty("line.separator"));

        final List<ACNVModeledSegment> segments = allSegments.stream()
                .filter(s -> !s.getContig().equals("2"))
                .filter(s -> !s.getContig().equals("6"))
                .filter(s -> s.getInterval().size() > lengthThreshold)
//                .filter(s -> s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_10) < log2crCredibleIntervalThreshold)
                .filter(s -> Double.isNaN(s.getMinorAlleleFractionPosteriorSummary().getCenter()) || s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_90) - s.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_10) < mafCredibleIntervalThreshold)
                        .collect(Collectors.toList());

        allSegments.stream().filter(s -> !segments.contains(s))
                .forEach(s -> {
                    output(writer, logger, "#fitered segment: " + s.getInterval() + " " +
                            s.getSegmentMeanPosteriorSummary().getCenter() + " " + s.getSegmentMeanPosteriorSummary().getLower() + " " + s.getSegmentMeanPosteriorSummary().getUpper() + " " +
                            s.getMinorAlleleFractionPosteriorSummary().getCenter() + " " + s.getMinorAlleleFractionPosteriorSummary().getLower() + " " + s.getMinorAlleleFractionPosteriorSummary().getUpper());
                    output(writer, logger, System.getProperty("line.separator"));
                });

        return segments;
    }

    private void outputModeller(final JavaSparkContext ctx,
                                final List<ACNVModeledSegment> segments,
                                final PloidyStatePrior variantPloidyStatePrior,
                                final int numPopulations,
                                final int numCells,
                                final int numSamples,
                                final int numBurnIn,
                                final TumorHeterogeneityModeller modeller,
                                final FileWriter writer,
                                final Logger logger) {
        //check statistics of global-parameter posterior samples (i.e., posterior mode and standard deviation)
        final Map<TumorHeterogeneityParameter, PosteriorSummary> globalParameterPosteriorSummaries =
                modeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);

        final double[] ploidySamples = Doubles.toArray(modeller.getPloidySamples());
        final double ploidyPosteriorMean = new Mean().evaluate(ploidySamples);
        final double ploidyPosteriorStandardDeviation = new StandardDeviation().evaluate(ploidySamples);
        output(writer, logger, "#ploidy: " + ploidyPosteriorMean + " " + ploidyPosteriorStandardDeviation);
        output(writer, logger, System.getProperty("line.separator"));

        final PosteriorSummary concentrationPosteriorSummary = globalParameterPosteriorSummaries.get(TumorHeterogeneityParameter.CONCENTRATION);
        final double concentrationPosteriorCenter = concentrationPosteriorSummary.getCenter();
        final double concentrationPosteriorStandardDeviation = (concentrationPosteriorSummary.getUpper() - concentrationPosteriorSummary.getLower()) / 2;
        output(writer, logger, "#concentration: " + concentrationPosteriorCenter + " " + concentrationPosteriorStandardDeviation);
        output(writer, logger, System.getProperty("line.separator"));

        final List<TumorHeterogeneityState.PopulationFractions> populationFractionsSamples = modeller.getPopulationFractionsSamples();
        final List<TumorHeterogeneityState.VariantProfileCollection> variantProfileCollectionSamples = modeller.getVariantProfileCollectionSamples();
        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            final int pi = populationIndex;
            final double[] populationFractionSamples = populationFractionsSamples.stream().mapToDouble(s -> s.get(pi)).toArray();
            final double populationFractionPosteriorMean = new Mean().evaluate(populationFractionSamples);
            final double populationFractionPosteriorStandardDeviation = new StandardDeviation().evaluate(populationFractionSamples);

            if (populationIndex != numPopulations - 1) {
                final double variantSegmentFraction = (double) IntStream.range(0, segments.size()).filter(i ->
                        (double) variantProfileCollectionSamples.stream().filter(vpc -> vpc.get(pi).isVariant(i)).count() / (numSamples - numBurnIn) > 0.9)
                        .count() / segments.size();
                output(writer, logger, "#population fraction " + populationIndex + ": " + populationFractionPosteriorMean + " " + populationFractionPosteriorStandardDeviation + " " + variantSegmentFraction);
                output(writer, logger, System.getProperty("line.separator"));
            } else {
                output(writer, logger, "#population fraction " + populationIndex + ": " + populationFractionPosteriorMean + " " + populationFractionPosteriorStandardDeviation);
                output(writer, logger, System.getProperty("line.separator"));
            }
        }

//        final List<TumorHeterogeneityState.PopulationIndicators> populationIndicatorsSamples = modeller.getPopulationIndicatorsSamples();
//        for (int cellIndex = 0; cellIndex < numCells; cellIndex++) {
//            final int ci = cellIndex;
//            final MultiSet<Integer> populationCounts = new HashMultiSet<>(populationIndicatorsSamples.stream().map(s -> s.get(ci)).collect(Collectors.toList()));
//            System.out.println("cell " + cellIndex + ": " + populationCounts);
//        }
//        System.out.println();

        //headers
        output(writer, logger, "POPULATION_INDEX\tSEGMENT_INDEX\tSEGMENT_INTERVAL\tIS_VARIANT_PROB\t");
        for (int variantPloidyStateIndex = 0; variantPloidyStateIndex < variantPloidyStatePrior.numPloidyStates(); variantPloidyStateIndex++) {
            final PloidyState variantPloidyState = variantPloidyStatePrior.ploidyStates().get(variantPloidyStateIndex);
            output(writer, logger, variantPloidyState.m() + "-" + variantPloidyState.n() + "\t");
        }
        output(writer, logger, System.getProperty("line.separator"));
        
        for (int populationIndex = 0; populationIndex < numPopulations; populationIndex++) {
            final int pi = populationIndex;
            final double[] populationFractionSamples = populationFractionsSamples.stream().mapToDouble(s -> s.get(pi)).toArray();
            final double populationFractionPosteriorMean = new Mean().evaluate(populationFractionSamples);

            if (populationIndex != numPopulations - 1 && populationFractionPosteriorMean >= 0.01) {
                for (int segmentIndex = 0; segmentIndex < segments.size(); segmentIndex++) {
                    final int si = segmentIndex;
                    final double isVariantPosteriorProbability = (double) variantProfileCollectionSamples.stream()
                            .filter(vpc -> vpc.get(pi).isVariant(si))
                            .count() / (numSamples - numBurnIn);
                    output(writer, logger, populationIndex + "\t" + segmentIndex + "\t" + segments.get(segmentIndex).getInterval() + "\t" + isVariantPosteriorProbability + "\t");

                    for (int variantPloidyStateIndex = 0; variantPloidyStateIndex < variantPloidyStatePrior.numPloidyStates(); variantPloidyStateIndex++) {
                        final int vpsi = variantPloidyStateIndex;
                        final double variantPloidyStateProbability = (double) variantProfileCollectionSamples.stream()
                                .filter(vpc -> vpc.get(pi).variantPloidyStateIndex(si) == vpsi)
                                .count() / (numSamples - numBurnIn);
                        output(writer, logger, String.format("%.3f", variantPloidyStateProbability));
                        if (variantPloidyStateIndex != variantPloidyStatePrior.numPloidyStates() - 1) {
                            output(writer, logger, "\t");
                        }
                    }
                    if (segmentIndex != segments.size() - 1) {
                        output(writer, logger, System.getProperty("line.separator"));
                    }
                }
            }
        }
    }
}