package analysis;

import htsjdk.samtools.SAMReadGroupRecord;

import java.io.IOException;
import java.io.File;

public interface AnonymizerAlgorithm {

    public final static String SHORT_READ_ALGORITHM = SAMReadGroupRecord.PlatformValue.ILLUMINA.name();

    public void anonymizeReads();

    public void queryReadsToExclude();

    public void setPartitions() throws IOException;
    public void setMinimumMappingQuality(int minimumMappingQuality);
    public void setThreadNumber(int threads);
    public void setMaxReadsInRam(int maxReadsInMemory);
    public void setTmpDir(File tmpDir);
    public void setMaxDepth(int maxDepth);
}
