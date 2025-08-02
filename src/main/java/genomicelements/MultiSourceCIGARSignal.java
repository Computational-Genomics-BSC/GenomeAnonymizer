package genomicelements;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class MultiSourceCIGARSignal implements ReadSignal{

    private Map<String, Integer> sources;
    private ReadSignal signal;
    private Signal.IndelSignalType indelType;

    public MultiSourceCIGARSignal(ReadSignal signal) {
        this.signal = signal;
        this.indelType = signal.getIndelSignalType();
        sources = new HashMap<>();
        this.addSource(signal.getReadAlnName(), signal.getInReadPosition());
    }

    public void addSource(String source, int inReadPosition) {
        sources.put(source, inReadPosition);
    }

    public Map<String, Integer> getSources() {
        return sources;
    }

    public Signal getReadAlnSignal(String readAlnName, int inReadPosition) {
        Signal answer;
        try {
            answer = (Signal) ((Signal) signal).clone();
            answer.setReadAlnName(readAlnName);
            answer.setInReadPosition(inReadPosition);
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }
        return answer;
    }

    public Signal.IndelSignalType getIndelType() {
        return indelType;
    }

    public void setIndelType(Signal.IndelSignalType indelType) {
        this.indelType = indelType;
    }

    public String getSequenceName() {
        return signal.getSequenceName();
    }

    public int getLocation() {
        return signal.getLocation();
    }

    public int getStart() {
        return signal.getStart();
    }

    public int getEnd() {
        return signal.getEnd();
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {
        signal.setSequenceIdx(sequenceIdx);
    }

    public int getSequenceIdx() {
        return signal.getSequenceIdx();
    }

    public int getLength() {
        return signal.getLength();
    }

    public String getReadAlnName() {
        return signal.getReadAlnName();
    }

    public int getInReadPosition() {
        return signal.getInReadPosition();
    }

    public boolean isFromNormalDataset() {
        return signal.isFromNormalDataset();
    }

    public boolean isGermline() {
        return signal.isGermline();
    }

    @Override
    public boolean isHandled() {
        return signal.isHandled();
    }

    @Override
    public void setHandled(boolean handled) {
        signal.setHandled(handled);
    }

    public boolean isClassifiedByDistance() {
        return signal.isClassifiedByDistance();
    }

    public byte[] getAltAllele() {
        return signal.getAltAllele();
    }

    @Override
    public boolean isSNV() {
        return false;
    }

    public boolean isIndel() {
        return signal.isIndel();
    }

    public boolean isSoftClip() {
        return signal.isSoftClip();
    }

    public boolean isDestructiveSignal() {
        return false;
    }

    public boolean isFromTumoralDataset() {
        return signal.isFromTumoralDataset();
    }

    public void setClassifiedByDistance(boolean classifiedByDistance) {
        signal.setClassifiedByDistance(classifiedByDistance);
    }

    public void setIsGermline(boolean isGermline) {
        signal.setIsGermline(isGermline);
    }

    @Override
    public void setIsFromNormalDataset(boolean isFromNormalDataset) {
        signal.setIsFromNormalDataset(isFromNormalDataset);
    }

    public String getSignalKey() {
        return signal.getSignalKey();
    }

    @Override
    public Signal.IndelSignalType getIndelSignalType() {
        return signal.getIndelSignalType();
    }

    @Override
    public Signal.Source getSource() {
        return signal.getSource();
    }

    @Override
    public boolean isIntraAlignmentSignal() {
        return true;
    }

    @Override
    public boolean isInterAlignmentSignal() {
        return false;
    }

    @Override
    public int compareTo(GenomicRegion genomicRegion) {
        return signal.compareTo(genomicRegion);
    }
}
