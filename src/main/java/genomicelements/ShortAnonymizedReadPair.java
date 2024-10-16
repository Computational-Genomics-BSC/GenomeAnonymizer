package genomicelements;

import htsjdk.samtools.fastq.FastqRecord;

public class ShortAnonymizedReadPair implements AnonymizedReadContainer{

    public static final int PAIR_1_IDX = 0;
    public static final int PAIR_2_IDX = 1;
    public static final String READ_PAIR_NAME_SEPARATOR = ";";

    private String readName;
    ShortAnonymizedRead pair1 = null;
    ShortAnonymizedRead pair2 = null;

    public ShortAnonymizedReadPair(ShortAnonymizedRead pair, boolean updatePairSequenceIfPossible) {
        readName = pair.getReadName();
        addOrUpdatePair(pair, updatePairSequenceIfPossible);
    }

    public ShortAnonymizedReadPair(ShortAnonymizedRead pair) {
        readName = pair.getReadName();
        addOrUpdatePair(pair);
    }

    @Override
    public boolean isWriteable() {
        if (pair1.isSupplementaryOrSecondary() || pair2.isSupplementaryOrSecondary()) return false;
        if (!pair1.isAnonymized() || !pair2.isAnonymized()) return false;
        return true;
    }

    @Override
    public FastqRecord[] getFastqRecords() {
        FastqRecord[] fqRecords = new FastqRecord[2];
        fqRecords[PAIR_1_IDX] = pair1.getFastqRecord();
        fqRecords[PAIR_2_IDX] = pair2.getFastqRecord();
        return fqRecords;
    }

    /**
     * Add or update pair from a new-found pair. This version is to be called with updatePairSequenceIfPossible = true,
     * for the anonymization phase
     * @param pair
     * @param updatePairSequenceIfPossible
     * @return true if the pair was added or updated, false if not
     */
    public boolean addOrUpdatePair(ShortAnonymizedRead pair, boolean updatePairSequenceIfPossible){
        int pairIdx = pair.getPairIdx();
        if(!hasPair(pairIdx)){
            addPair(pair);
            return true;
        }
        else{
            if (updatePairSequenceIfPossible){
                ShortAnonymizedRead savedPair = getPair(pairIdx);
                if (savedPair.isSupplementaryOrSecondary() && !pair.isSupplementaryOrSecondary()){
                    savedPair.setSequenceArray(pair.getSequenceArray(), pair.getQualitiesArray());
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Add or update pair from a new-found pair. This version is to be called for the variation discovery phase
     * @param pair
     * @return true if the pair was added or updated, false if not
     */
    public boolean addOrUpdatePair(ShortAnonymizedRead pair){
        return addOrUpdatePair(pair, false);
    }

    public boolean hasBothPairs(){
        return (pair1 != null && pair2 != null);
    }

    public boolean hasPair(int pairIdx){
        if (pairIdx == PAIR_1_IDX){
            return hasPair1();
        }
        if (pairIdx == PAIR_2_IDX){
            return hasPair2();
        }
        return false;
    }

    public void addPair(ShortAnonymizedRead pair){
        int pairIdx = pair.getPairIdx();
        if (pairIdx == PAIR_1_IDX){
            addPair1(pair);
        }
        if (pairIdx == PAIR_2_IDX){
            addPair2(pair);
        }
    }

    public void addPair1(ShortAnonymizedRead pair1){
        this.pair1 = pair1;
    }

    public void addPair2(ShortAnonymizedRead pair2){
        this.pair2 = pair2;
    }

    public boolean hasPair1(){
        return pair1 != null;
    }

    public boolean hasPair2(){
        return pair2 != null;
    }

    public String getReadName() {
        return readName;
    }

    public ShortAnonymizedRead getPair(int pairIdx){
        ShortAnonymizedRead answer = null;
        if (pairIdx == PAIR_1_IDX){
            answer = pair1;
        }
        else if (pairIdx == PAIR_2_IDX){
            answer =  pair2;
        }
        return answer;
    }

    public ShortAnonymizedRead getPair1() {
        return pair1;
    }

    public ShortAnonymizedRead getPair2() {
        return pair2;
    }


}
