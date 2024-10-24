package analysis;

import genomicelements.CalledVariation;
import genomicelements.ShortAnonymizedReadPair;
import htsjdk.samtools.SAMReadGroupRecord;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public interface AnonymizerAlgorithm {

    public final static String SHORT_READ_ALGORITHM = SAMReadGroupRecord.PlatformValue.ILLUMINA.name();

    //public Map<String, Map<Integer, List<CalledVariation>>> callVariation(String normalPath, String tumorPath, String refGenome, String mode, String vcfFile)throws IOException;

    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed)throws IOException;

    public Map<String, ShortAnonymizedReadPair> getAnonymizedReadContainer();

    public void writeUnmodifiedReads(String normalPath, String tumorPath, String outputPrefix, boolean compressed)throws IOException;

    public void setRemoveUnmapped(boolean removeUnmapped);
}
