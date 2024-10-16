package genomicelements;

import htsjdk.samtools.fastq.FastqRecord;

public interface AnonymizedReadContainer {

    public String getReadName();

    public boolean isWriteable();

    public FastqRecord[] getFastqRecords();

}
