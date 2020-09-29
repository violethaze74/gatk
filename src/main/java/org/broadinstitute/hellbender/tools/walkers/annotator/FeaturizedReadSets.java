package org.broadinstitute.hellbender.tools.walkers.annotator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * For each sample and for each allele a list feature vectors of supporting reads
 * In order to reduce the number of delimiter characters, we flatten featurized reads.  For example, suppose allele 1 has featurized reads
 * [1,2] and [3,4] and allele 2 has featurized reads [5,6] and [7,8], the annotation is
 * 1,2,3,4|5,6,7,8
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY,
        summary="Featurized read sets for Mutect3 training data")
public class FeaturizedReadSets extends GenotypeAnnotation {
    public static final int DEFAULT_BASE_QUALITY = 25;

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        Utils.nonNull(g);
        Utils.nonNull(gb);

        if ( likelihoods == null) {
            return;
        }

        final Map<Allele, List<Integer>> flattenedFeaturizedReads = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        Utils.stream(likelihoods.bestAllelesBreakingTies())
                .filter(ba -> ba.isInformative())
                .forEach(ba -> flattenedFeaturizedReads.get(ba.allele).addAll(featurize(ba.evidence, vc)));

        final List<List<Integer>> dataInAlleleOrder = likelihoods.alleles().stream().map(flattenedFeaturizedReads::get).collect(Collectors.toList());
        final List<String> stringsInAlleleOrder = dataInAlleleOrder.stream().map(list -> StringUtils.join(list, ",")).collect(Collectors.toList());

        final String annotation = AnnotationUtils.encodeAnyASListWithRawDelim(stringsInAlleleOrder);


        gb.attribute(GATKVCFConstants.FEATURIZED_READ_SETS_KEY, annotation);
    }



    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.FEATURIZED_READ_SETS_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getFormatLine(getKeyNames().get(0)));
    }

    private List<Integer> featurize(final GATKRead read, final VariantContext vc) {
        final List<Integer> result = new ArrayList<>();
        result.add(read.getMappingQuality());
        result.add(BaseQuality.getBaseQuality(read, vc).orElse(DEFAULT_BASE_QUALITY));
        result.add(read.isFirstOfPair() ? 1 : 0);
        result.add(read.isReverseStrand() ? 1 : 0);
        result.add(ReadPosition.getPosition(read, vc).orElse(0));
        result.add(Math.abs(read.getFragmentLength()));

        return result;
    }

}
