package analysis;

import genomicelements.CalledVariation;
import htsjdk.tribble.SimpleFeature;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
    private SimpleFeature region;

    public MultithreadClassifier(String normalPath, String tumorPath, String refGenome, String mode, String vcfFile, SimpleFeature region){
        this.normalPath = normalPath;
        this.tumorPath = tumorPath;
        this.refGenome = refGenome;
        this.mode = mode;
        this.vcfFile = vcfFile;
        this.region = region;
    }

    @Override
    public void run() {
        try {
            classifier.callVariation(normalPath, tumorPath, refGenome, mode, vcfFile, region);
        } catch (Exception e) {
            LOGGER.severe("Exception in thread for region: " + region.getContig() + " " + region.getStart() + " " + region.getEnd() +
                    " halting execution prematurely");
            LOGGER.severe(e.getMessage());
            //System.exit(1);
        }
    }

    public String getContig(){
        return this.region.getContig();
    }

    public Map<String, Map<Integer, List<CalledVariation>>> getAnswer(){
        assert(classifier.getPotentialGermlinesPerRead().size() == 1): "The result of this classifier is incorrect: "
                + region.getContig() + " " + region.getStart() + " " + region.getEnd();
        //return classifier.getPotentialGermlinesPerRead().getOrDefault(region.getContig(), new HashMap<>());
        return classifier.getPotentialGermlinesPerRead();
    }
}
