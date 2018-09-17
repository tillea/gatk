package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Given annotated intervals output by {@link AnnotateIntervals} and/or counts collected on those intervals output
 * by {@link CollectReadCounts}, outputs a filtered interval list.  Parameters for filtering based on the annotations
 * and counts can be adjusted.  Annotation-based filters will be applied first, followed by count-based
 * filters.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
           Intervals to be filtered.
 *         The argument {@code interval-merging-rule} must be set to {@link IntervalMergingRule#OVERLAPPING_ONLY}
 *         and all other common arguments for interval padding or merging must be set to their defaults.
 *         A blacklist of regions in which intervals should always be filtered (regardless of other annotation-based
 *         or count-based filters) may also be provided via -XL; this can be used to filter pseudoautosomal regions
 *         (PARs), for example.
 *     </li>
 *     <li>
 *         (Optional) Annotated-intervals file from {@link AnnotateIntervals}.
 *         Must contain the intervals to be filtered as a subset.  Must be provided if no counts files are provided.
 *     </li>
 *     <li>
 *         (Optional) Counts files (TSV or HDF5 output of {@link CollectReadCounts}).
 *         Must contain the intervals to be filtered as a subset.  Must be provided if no annotated-intervals file
 *         is provided.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Filtered Picard interval-list file.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk FilterIntervals \
 *          -L intervals.interval_list \
 *          -XL blacklist_intervals.interval_list \
 *          -I sample_1.counts.hdf5 \
 *          -I sample_2.counts.hdf5 \
 *          ... \
 *          --annotated-intervals annotated_intervals.tsv \
 *          -O filtered_intervals.interval_list
 * </pre>
 *
 * <pre>
 *     gatk FilterIntervals \
 *          -L intervals.interval_list \
 *          --annotated-intervals annotated_intervals.tsv \
 *          -O filtered_intervals.interval_list
 * </pre>
 *
 * <pre>
 *     gatk FilterIntervals \
 *          -L intervals.interval_list \
 *          -I sample_1.counts.hdf5 \
 *          -I sample_2.counts.hdf5 \
 *          ... \
 *          -O filtered_intervals.interval_list
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Filters intervals based on annotations and/or count statistics",
        oneLineSummary = "Filters intervals based on annotations and/or count statistics",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class FilterIntervals extends CommandLineProgram {
    public static final String MINIMUM_GC_CONTENT_LONG_NAME = "minimum-gc-content";
    public static final String MAXIMUM_GC_CONTENT_LONG_NAME = "maximum-gc-content";
    public static final String MINIMUM_MAPPABILITY_LONG_NAME = "minimum-mappability";
    public static final String MAXIMUM_MAPPABILITY_LONG_NAME = "maximum-mappability";
    public static final String MINIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME = "minimum-segmental-duplication-content";
    public static final String MAXIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME = "maximum-segmental-duplication-content";

    @Argument(
            doc = "Input file containing annotations for genomic intervals (output of AnnotateIntervals).  " +
                    "All intervals specified via -L must be contained.  " +
                    "Must be provided if no counts files are provided.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAnnotatedIntervalsFile = null;

    @Argument(
            doc = "Input TSV or HDF5 files containing integer read counts in genomic intervals (output of CollectReadCounts).  " +
                    "All intervals specified via -L must be contained.  " +
                    "Must be provided if no annotated-intervals file is provided.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Output file for filtered intervals.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFilteredIntervalsFile;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection
            = new OptionalIntervalArgumentCollection();

    @Argument(
            doc = "Minimum allowed value for GC-content annotation (inclusive).",
            fullName = MINIMUM_GC_CONTENT_LONG_NAME,
            optional = true
    )
    private double minimumGCContent = 0.1;

    @Argument(
            doc = "Maximum allowed value for GC-content annotation (inclusive).",
            fullName = MAXIMUM_GC_CONTENT_LONG_NAME,
            optional = true
    )
    private double maximumGCContent = 0.9;

    @Argument(
            doc = "Minimum allowed value for mappability annotation (inclusive).",
            fullName = MINIMUM_MAPPABILITY_LONG_NAME,
            optional = true
    )
    private double minimumMappability = 0.9;

    @Argument(
            doc = "Maximum allowed value for mappability annotation (inclusive).",
            fullName = MAXIMUM_MAPPABILITY_LONG_NAME,
            optional = true
    )
    private double maximumMappability = 1.;

    @Argument(
            doc = "Minimum allowed value for segmental-duplication-content annotation (inclusive).",
            fullName = MINIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME,
            optional = true
    )
    private double minimumSegmentalDuplicationContent = 0.;

    @Argument(
            doc = "Maximum allowed value for segmental-duplication-content annotation (inclusive).",
            fullName = MAXIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME,
            optional = true
    )
    private double maximumSegmentalDuplicationContent = 0.5;

    private SimpleIntervalCollection specifiedIntervals;
    private AnnotatedIntervalCollection annotatedIntervals;
    private RealMatrix readCountMatrix;

    @Override
    public Object doWork() {
        validateArguments();

        if (inputReadCountFiles.isEmpty()) {
            //only annotated intervals provided (no counts)
            //we end up reading the annotated intervals twice, first to get the sequence dictionary; probably not worth cleaning up
            annotatedIntervals = new AnnotatedIntervalCollection(inputAnnotatedIntervalsFile);
            specifiedIntervals = new SimpleIntervalCollection(
                    annotatedIntervals.getMetadata(),
                    intervalArgumentCollection.getIntervals(annotatedIntervals.getMetadata().getSequenceDictionary()));
            Utils.validateArg(specifiedIntervals.size() != 0, "At least one interval must be specified.");
            Utils.validateArg(new HashSet<>(annotatedIntervals.getIntervals()).containsAll(specifiedIntervals.getIntervals()),
                    "Annotated intervals do not contain all specified intervals.");
        } else {
            //counts provided
            final File firstReadCountFile = inputReadCountFiles.get(0);
            specifiedIntervals = CopyNumberArgumentValidationUtils.resolveIntervals(
                    firstReadCountFile, intervalArgumentCollection, logger);
            if (inputAnnotatedIntervalsFile != null) {
                //both annotated intervals and counts provided
                annotatedIntervals = CopyNumberArgumentValidationUtils.validateAnnotatedIntervalsSubset(
                        inputAnnotatedIntervalsFile, specifiedIntervals, logger);
            }
            readCountMatrix = constructReadCountMatrix(inputReadCountFiles);
        }
        final SimpleIntervalCollection filteredIntervals = filterIntervals();
        logger.info(String.format("Writing filtered intervals to %s...", outputFilteredIntervalsFile));
        filteredIntervals.write(outputFilteredIntervalsFile);

        return "SUCCESS";
    }

    private void validateArguments() {
        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);
        if (inputAnnotatedIntervalsFile == null && inputReadCountFiles.isEmpty()) {
            throw new UserException("Must provide annotated intervals or counts.");
        }
        if (inputAnnotatedIntervalsFile != null) {
            IOUtils.canReadFile(inputAnnotatedIntervalsFile);
        }
        inputReadCountFiles.forEach(IOUtils::canReadFile);
        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-count files cannot contain duplicates.");
    }

    private RealMatrix constructReadCountMatrix(final List<File> inputReadCountFiles) {
        logger.info("Validating and aggregating input read-counts files...");
        final int numSamples = inputReadCountFiles.size();
        final int numIntervals = specifiedIntervals.size();
        final Set<SimpleInterval> intervalSubset = new HashSet<>(specifiedIntervals.getRecords());
        final RealMatrix readCountMatrix = new Array2DRowRealMatrix(numSamples, numIntervals);
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-counts file %s (%d / %d)", inputReadCountFile, sampleIndex + 1, numSamples));
            final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);
            if (!CopyNumberArgumentValidationUtils.isSameDictionary(
                    readCounts.getMetadata().getSequenceDictionary(),
                    specifiedIntervals.getMetadata().getSequenceDictionary())) {
                logger.warn(String.format("Sequence dictionary for read-counts file %s does not match those in other read-counts files.", inputReadCountFile));
            }
            final double[] subsetReadCounts = readCounts.getRecords().stream()
                    .filter(c -> intervalSubset.contains(c.getInterval()))
                    .mapToDouble(SimpleCount::getCount)
                    .toArray();
            Utils.validateArg(subsetReadCounts.length == intervalSubset.size(),
                    String.format("Intervals for read-count file %s do not contain all specified intervals.",
                            inputReadCountFile));
            readCountMatrix.setRow(sampleIndex, subsetReadCounts);
        }
        return readCountMatrix;
    }

    private SimpleIntervalCollection filterIntervals() {
        //filter by annotations
        List<AnnotatedInterval> filteredAnnotatedIntervals = annotatedIntervals.getRecords();
        final List<AnnotationKey<?>> annotationKeys = annotatedIntervals.getRecords().get(0).getAnnotationMap().getKeys();
        if (annotationKeys.contains(CopyNumberAnnotations.GC_CONTENT)) {    //this should always be true, but we check it anyway
            filteredAnnotatedIntervals = filteredAnnotatedIntervals.stream()
                    .filter(i -> {
                        final double value = i.getAnnotationMap().getValue(CopyNumberAnnotations.GC_CONTENT);
                        return minimumGCContent <= value && value <= maximumGCContent;})
                    .collect(Collectors.toList());
        }
        if (annotationKeys.contains(CopyNumberAnnotations.MAPPABILITY)) {
            filteredAnnotatedIntervals = filteredAnnotatedIntervals.stream()
                    .filter(i -> {
                        final double value = i.getAnnotationMap().getValue(CopyNumberAnnotations.MAPPABILITY);
                        return minimumMappability <= value && value <= maximumMappability;})
                    .collect(Collectors.toList());
        }
        if (annotationKeys.contains(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT)) {
            filteredAnnotatedIntervals = filteredAnnotatedIntervals.stream()
                    .filter(i -> {
                        final double value = i.getAnnotationMap().getValue(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT);
                        return minimumSegmentalDuplicationContent <= value && value <= maximumSegmentalDuplicationContent;})
                    .collect(Collectors.toList());
        }
        return new SimpleIntervalCollection(
                specifiedIntervals.getMetadata(),
                filteredAnnotatedIntervals.stream().map(AnnotatedInterval::getInterval).collect(Collectors.toList()));
    }
}
