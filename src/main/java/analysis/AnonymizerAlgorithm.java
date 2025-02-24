package analysis;

import genomicelements.GenomicRegion;
import genomicelements.Signal;
import htsjdk.samtools.SAMReadGroupRecord;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

public interface AnonymizerAlgorithm {

    public final static String SHORT_READ_ALGORITHM = SAMReadGroupRecord.PlatformValue.ILLUMINA.name();

    public void queryReadsToExclude(String normalPath, String tumorPath, int threads) throws Exception;
    public Set<String> getReadsToExclude();
    public void setGenomicPartitions(List<GenomicRegion> genomicPartitions);
    public void setReadGermlinesToAnonymize(Map<String, List<Signal>> readGermlinesToAnonymize);
    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed)throws IOException;

    public void setQueryRegions(List<GenomicRegion> regions);
    public void setCanvasFiles(String normalCanvasFileName, String tumoralCanvasFileName);
}
