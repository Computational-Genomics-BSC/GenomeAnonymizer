package genomicelements;

public class Signature implements GenomicRegion{

    private String sequenceName;
    private int start;
    private int end;
    private int sequenceIdx;
    private int length;
    private byte[] sequenceBytes;
    private String source;


    @Override
    public String getSequenceName() {
        return "";
    }

    @Override
    public int getSequenceIdx() {
        return 0;
    }

    @Override
    public int getStart() {
        return 0;
    }

    @Override
    public int getEnd() {
        return 0;
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {

    }

    public enum Source{
        SOFT_CLIP(1, "SOFT_CLIP"),
        INSERT_SIZE(2, "INSERT_SIZE");

        private final int value;
        private final String name;

        Source(int value, String name){
            this.name = name;
            this.value = value;
        }

        public int getValue() {
            return value;
        }

        public String getName() {
            return name;
        }
    }
}
