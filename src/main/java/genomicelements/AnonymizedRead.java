package genomicelements;

import htsjdk.samtools.SAMRecord;

import java.util.List;
import java.util.Map;

// May be better as a parent class
public interface AnonymizedRead{

    public String getReadAlignmentId();

    public void setVariantsToAnonymize(List<CalledVariation> variants);

    public void anonymizeVariants();

    public boolean isAnonymized();

    public SAMRecord getAnonymizedSamRecord();
}
