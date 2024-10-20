package analysis;

import genomicelements.ShortAnonymizedReadPair;
import htsjdk.samtools.SAMReadGroupRecord;

import java.io.IOException;
import java.util.Map;

public interface AnonymizerAlgorithm {

    public final static String SHORT_READ_ALGORITHM = SAMReadGroupRecord.PlatformValue.ILLUMINA.name();

    public void callVariation(String normalPath, String tumorPath, String refGenome)throws IOException;

    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed)throws IOException;

    public Map<String, ShortAnonymizedReadPair> getAnonymizedReadContainer();

    public void writeUnmodifiedReads(String normalPath, String tumorPath, String outputPrefix, boolean compressed)throws IOException;

    public void setRemoveUnmapped(boolean removeUnmapped);
}
