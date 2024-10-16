package genomicelements;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Locatable;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

// May be better as a parent class
public interface AnonymizedRead extends Locatable {

    public byte[] getSequenceArray();

    public Map<String, List<CalledVariation>> getVariantsToAnonymize();

    public void modifyBaseInRead(int inReadPosition, byte asciiBase);

    public boolean addVariantToAnonymize(CalledVariation variation);

    public void anonymizeVariantsInRead();

    public String getUniqueReadName();

    public boolean isAnonymized();

    public FastqRecord getFastqRecord();
}
