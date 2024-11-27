package genomicelements;

import htsjdk.samtools.SAMRecord;

import java.util.List;

public interface AnonymizedRead{

    public String getReadAlignmentId();

    public void setVariantsToAnonymize(List<CalledVariation> variants);

    public void anonymizeVariants();

    public boolean isAnonymized();

    public SAMRecord getAnonymizedSamRecord();
}
