package genomicelements;

import htsjdk.samtools.SAMRecord;

import java.util.List;

public interface AnonymizedRead {
    public String getSequenceName();
    public int getStart();
    public int getEnd();
    public boolean addSignalToAnonymize(Signal signal);
    public String getReadName();
    public String getReadAlignmentId();
    public void setSignalsToAnonymize(List<Signal> signals);
    public void anonymizeRead();
    public boolean isFromNormalDataset();
    public boolean isFromTumoralDataset();
    public boolean isAnonymized();
    public SAMRecord getAnonymizedSamRecord();
    public void setReferenceContigSequence(byte[] referenceContigSequence);
    public void rescueIndelSignalToAnonymize(Signal signal);
}
