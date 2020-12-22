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
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
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

    private static final int DEFAULT_MAX_REF_COUNT = Integer.MAX_VALUE;

    private static final int FEATURES_PER_READ = 9;

    // downsample ref reads to this count if needed
    private final int maxRefCount;

    public FeaturizedReadSets(final int maxRefCount) {
        this.maxRefCount = maxRefCount;
    }

    public FeaturizedReadSets() {
        this(DEFAULT_MAX_REF_COUNT);
    }

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



        final Map<Allele, List<GATKRead>> readsByAllele = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        Utils.stream(likelihoods.bestAllelesBreakingTies())
                .filter(ba -> ba.isInformative())
                .forEach(ba -> readsByAllele.get(ba.allele).add(ba.evidence));

        // downsample if necessary
        final Allele refAllele = likelihoods.alleles().stream().filter(Allele::isReference).findFirst().get();
        if (readsByAllele.get(refAllele).size() > maxRefCount) {
            Collections.shuffle(readsByAllele.get(refAllele));
            readsByAllele.put(refAllele, readsByAllele.get(refAllele).subList(0, maxRefCount));
        }

        final List<String> stringsInAlleleOrder = vc.getAlleles().stream()
                .map(allele -> {
                            final List<GATKRead> reads = readsByAllele.get(allele);
                            final List<Integer> flattened = new ArrayList<>(reads.size() * FEATURES_PER_READ);
                            reads.forEach(read -> flattened.addAll(featurize(read, vc)));
                            return StringUtils.join(flattened, ",");
                        }).collect(Collectors.toList());


        final String annotation = AnnotationUtils.encodeAnyASListWithRawDelim(stringsInAlleleOrder);


        gb.attribute(GATKVCFConstants.FEATURIZED_READ_SETS_KEY, annotation);
    }



    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.FEATURIZED_READ_SETS_KEY);
    }

    private List<Integer> featurize(final GATKRead read, final VariantContext vc) {
        final List<Integer> result = new ArrayList<>();
        result.add(read.getMappingQuality());
        result.add(BaseQuality.getBaseQuality(read, vc).orElse(DEFAULT_BASE_QUALITY));
        result.add(read.isFirstOfPair() ? 1 : 0);
        result.add(read.isReverseStrand() ? 1 : 0);

        // distances from ends of read
        final int readPosition = ReadPosition.getPosition(read, vc).orElse(0);
        result.add(readPosition);
        result.add(read.getLength() - readPosition);


        result.add(Math.abs(read.getFragmentLength()));

        // distances from ends of fragment
        final int fragmentStart = Math.min(read.getMateStart(), read.getUnclippedStart());
        final int fragmentEnd = fragmentStart + Math.abs(read.getFragmentLength());
        result.add(vc.getStart() - fragmentStart);
        result.add(fragmentEnd - vc.getEnd());

        // TODO: mismatches versus best haplotype



        Utils.validate(result.size() == FEATURES_PER_READ, "Wrong number of features");

        return result;
    }

}
