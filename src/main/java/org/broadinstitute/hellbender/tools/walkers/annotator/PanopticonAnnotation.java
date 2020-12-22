package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Annotations relevant to the INFO field of the variant file (ie annotations for sites).
 */
public abstract class PanopticonAnnotation extends VariantAnnotation{

    // must implement:
    // public abstract List<String> getKeyNames()
    // from VariantAnnotation superclass

    /**
     * Computes the annotation for the given genotype and the likelihoods per read.
     * Expected to modified the passed genotype builder.
     *
     * @param ref Reference context, may be null
     * @param vc Variant to be annotated. Not null.
     * @param likelihoods matrix of likelihoods indexed by allele and read
     * @param g the genotype to annotate. May be null.
     * @param gb the builder to modify and annotations to. Not null.
     */
    public abstract void annotate(final ReferenceContext ref,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final AlleleLikelihoods<GATKRead, Allele> likelihoods);

    /**
     * Computes the annotation for the given variant and the likelihoods per read.
     * Returns a map from annotation keys to values (may be empty if no annotation is to be added).
     *
     * @param ref Reference context, may be null
     * @param vc Variant to be annotated. Not null.
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     */
    public abstract Map<String, Object> annotate(final ReferenceContext ref,
                                                 final VariantContext vc,
                                                 final AlleleLikelihoods<GATKRead, Allele> likelihoods);
}