package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.util.Lazy;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Superclass of all variant annotations.
 */
public abstract class VariantAnnotation implements Annotation{

    // is this an INFO or FORMAT annotation
    public abstract VCFCompoundHeaderLine.SupportedHeaderLineType annotationType();

    private final Lazy<List<VCFCompoundHeaderLine>> headerLines = new Lazy<>(() -> getKeyNames().stream().map(key -> {
        switch (annotationType()) {
            case INFO:
                return GATKVCFHeaderLines.getInfoLine(key, true);
            case FORMAT:
                return GATKVCFHeaderLines.getFormatLine(key, true);
            default:
                throw new IllegalStateException("Unsupported annotation type: " + annotationType());

        }}).collect(Collectors.collectingAndThen(Collectors.toList(), ImmutableList::copyOf)));

    // Return the descriptions used for the VCF INFO or FORMAT meta field.
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return headerLines.get();
    }

    /**
     * Return the keys
     */
    public abstract List<String> getKeyNames();

    @Override
    public String toString() {
        return getClass().getSimpleName();
    }
}