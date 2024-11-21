package analysis;

import genomicelements.CalledVariation;
import genomicelements.ShortAnonymizedReadPair;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.tribble.SimpleFeature;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

public interface AnonymizerAlgorithm {

    public final static String SHORT_READ_ALGORITHM = SAMReadGroupRecord.PlatformValue.ILLUMINA.name();

    public void queryReadsToExclude(String normalPath, String tumorPath, List<SimpleFeature> partitions, int threads) throws Exception;

    public void setReadGermlinesToAnonymize(Map<String, Map<Integer, List<CalledVariation>>> readGermlinesToAnonymize);

    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed)throws IOException;

    public void writeUnmodifiedReads(String normalPath, String tumorPath, String outputPrefix, boolean compressed)throws IOException;

    public Set<String> getReadsToExclude();
}
