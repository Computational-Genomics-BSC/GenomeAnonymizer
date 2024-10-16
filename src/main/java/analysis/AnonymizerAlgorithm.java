package analysis;

import genomicelements.AnonymizedReadContainer;
import genomicelements.PairedPileup;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.IOException;
import java.util.Map;

public interface AnonymizerAlgorithm {

    public final static String SHORT_READ_ALGORITHM = SAMReadGroupRecord.PlatformValue.ILLUMINA.name();

    public void anonymizeReads(Iterable<PairedPileup> readFileIterable, IndexedFastaSequenceFile referenceFasta) throws IOException;

    public Map<String, AnonymizedReadContainer> getAnonymizedReadContainer();
}
