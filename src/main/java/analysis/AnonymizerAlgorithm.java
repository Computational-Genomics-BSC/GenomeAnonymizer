package analysis;

import genomicelements.GenomicRegion;
import genomicelements.PairCalledVariation;
import genomicelements.Signal;
import htsjdk.samtools.SAMReadGroupRecord;

import java.io.IOException;
import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

public interface AnonymizerAlgorithm {

    public final static String SHORT_READ_ALGORITHM = SAMReadGroupRecord.PlatformValue.ILLUMINA.name();

    public void anonymizeReads();

    public void queryReadsToExclude();

    public void setPartitions() throws IOException;
    public void setMinimumMappingQuality(int minimumMappingQuality);
    public void setThreadNumber(int threads);
    public void setMaxReadsInRam(int maxReadsInMemory);
    public void setTmpDir(File tmpDir);

    //TODO: Decide if this stays or not
    public void setVCFVariantsToKeep(Map<String, Map<Integer, PairCalledVariation>> variantsToKeep);
}
