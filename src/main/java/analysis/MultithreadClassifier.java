package analysis;

import genomicelements.CalledVariation;
import genomicelements.GenomicRegion;
import genomicelements.Signal;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Helper class to multithread the variation calling phase
 * @author Nicolas Gaitan
 */
public class MultithreadClassifier implements Runnable{

    private static final Logger LOGGER = Logger.getLogger(MultithreadClassifier.class.getName());

    private final VariationClassifier classifier = new VariationClassifier();
    private String normalPath;
    private String tumorPath;
    private String refGenome;
    private String mode;
    private String vcfFile;
    private GenomicRegion region;
    private Set<String> readsToExclude;

    public MultithreadClassifier(String normalPath, String tumorPath, String refGenome,
                                 GenomicRegion region, Set<String> readsToExclude){
        this.normalPath = normalPath;
        this.tumorPath = tumorPath;
        this.refGenome = refGenome;
        this.region = region;
        this.readsToExclude = readsToExclude;
    }

    @Override
    public void run() {
        try {
            if(!readsToExclude.isEmpty()) classifier.setReadsToExclude(readsToExclude);
            long startcallVariation = System.currentTimeMillis();
            classifier.callVariation(normalPath, tumorPath, refGenome, region);
            long endcallVariation = System.currentTimeMillis();
            classifier.METHOD_TIME_MAP.put("callVariation", endcallVariation-startcallVariation);
            LOGGER.info("Finished variation analysis of genomic region: SEQ=" + region.getSequenceName() + " POS=" + region.getStart() + " END=" + region.getEnd());
            //DEBUG
            StringBuilder msg = new StringBuilder();
            for (Map.Entry<String, Long> entry : classifier.METHOD_TIME_MAP.entrySet()){
                msg.append(entry.getKey()).append(": ");
                long timeInSeconds = entry.getValue();
                msg.append(timeInSeconds).append("\n");
            }
            LOGGER.info("PARTITION TIME TABLE: " + "\n" + msg);
            //DEBUG
        } catch (Exception e) {
            LOGGER.log(Level.SEVERE, "Exception in thread classifying variants in region: "
                            + region.getSequenceName() + " " + region.getStart() + " " + region.getEnd() +
                            " halting execution prematurely",
                    e);
            System.exit(1);
        }
    }

    public String getSequenceName(){
        return this.region.getSequenceName();
    }
    public int getStart(){return this.region.getStart();}
    public int getEnd(){return this.region.getEnd();}


    public Map<String, List<Signal>>  getAnswer(){
        assert(classifier.getPotentialGermlinesPerRead().size() == 1): "The result of this classifier is incorrect: "
                + region.getSequenceName() + " " + region.getStart() + " " + region.getEnd();
        return classifier.getPotentialGermlinesPerRead();
    }

    public void setVCFVariantsToKeep(Map<String, Map<Integer,CalledVariation>> somaticVariantsToKeep){
        classifier.setVCFVariantsToKeep(somaticVariantsToKeep);
    }
}
