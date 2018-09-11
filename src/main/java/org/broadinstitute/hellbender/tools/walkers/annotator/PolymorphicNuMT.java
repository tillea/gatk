package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of alts indicates it could be an autosomal false positive.")
public class PolymorphicNuMT extends GenotypeAnnotation implements Annotation {
    protected final OneShotLogger warning = new OneShotLogger(this.getClass());
    public static final String POTENTIAL_POLYMORPHIC_NUMT = "POTENTIAL_POLYMORPHIC_NUMT";
    private static final double LOWER_BOUND_PROB = .1;
    private int MIN_AUTOSOMAL_HET;
    private int MAX_AUTOSOMAL_HET;
    private int MIN_AUTOSOMAL_HOM_ALT;
    private int MAX_AUTOSOMAL_HOM_ALT;

    public PolymorphicNuMT(final double lambda){
        PoissonDistribution autosomalCoverage = new PoissonDistribution(lambda);
        this.MIN_AUTOSOMAL_HOM_ALT = autosomalCoverage.inverseCumulativeProbability(LOWER_BOUND_PROB);
        this.MAX_AUTOSOMAL_HOM_ALT = autosomalCoverage.inverseCumulativeProbability(1 - LOWER_BOUND_PROB);
        this.MIN_AUTOSOMAL_HET = MIN_AUTOSOMAL_HOM_ALT / 2;
        this.MAX_AUTOSOMAL_HET = MAX_AUTOSOMAL_HOM_ALT / 2;
    }

    // Barclay requires that each annotation define a constructor that takes no arguments
    public PolymorphicNuMT(){ }

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY, () -> null, -1);
        if (tumorLods==null) {
            warning.warn("One or more variant contexts is missing the 'TLOD' annotation, StrandArtifact will not be computed for these VariantContexts");
            return;
        }
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxTumorLod);
        Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());
        final long numAltReads = bestAlleles.stream().filter(ba -> ba.isInformative() && ba.allele.equals(altAlelle)).count();
        if ( (numAltReads > MIN_AUTOSOMAL_HET && numAltReads < MAX_AUTOSOMAL_HET) || (numAltReads > MIN_AUTOSOMAL_HOM_ALT && numAltReads < MAX_AUTOSOMAL_HOM_ALT) ) {
            gb.attribute(POTENTIAL_POLYMORPHIC_NUMT, "true");
        }
    }
    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(POTENTIAL_POLYMORPHIC_NUMT, 1, VCFHeaderLineType.String, "Potentially a polymorphic NuMT false positive rather than a real mitochondrial variant."));
    }
    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(POTENTIAL_POLYMORPHIC_NUMT);
    }
}
