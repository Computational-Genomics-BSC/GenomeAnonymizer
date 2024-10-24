package analysis;

import genomicelements.CalledVariation;
import htsjdk.samtools.util.IntervalList;
import htsjdk.tribble.SimpleFeature;

import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * Helper class to multithread the variation calling phase
 * @author Nicolas Gaitan
 */
public class MultithreadClassifier implements Runnable{

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
            throw new RuntimeException(e);
        }
    }

    public Map<String, Map<Integer, List<CalledVariation>>> getResult(){
        return classifier.getPotentialGermlinesPerRead();
    }
}
