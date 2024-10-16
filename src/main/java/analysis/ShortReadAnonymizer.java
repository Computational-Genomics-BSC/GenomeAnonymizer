package analysis;

import genomicelements.AnonymizedReadContainer;
import genomicelements.CalledVariation;
import genomicelements.PairedPileup;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.IOException;
import java.util.*;

public class ShortReadAnonymizer implements AnonymizerAlgorithm{

    public static final int SLIDING_WINDOW_LIMIT = 200;

    Iterable<PairedPileup> reader;
    Map<String, AnonymizedReadContainer> anonReadContainer;
    IndexedFastaSequenceFile referenceWalker;

    public ShortReadAnonymizer() {
        anonReadContainer = new HashMap<>();;
    }

    @Override
    public void anonymizeReads(Iterable<PairedPileup> pairPileupReader, IndexedFastaSequenceFile referenceWalker) throws IOException {
        Map<Integer, List<CalledVariation>> variationPerPos = new HashMap<>();
        Set<String> seenReads = new HashSet<>();
        VariationClassifier classifier = new VariationClassifier();
        int p = 1;
        for (PairedPileup pileup : pairPileupReader){
            int pos = pileup.getReferencePos();
            classifier.classifyVariationInPairedPileup(variationPerPos, pileup, seenReads, referenceWalker);
            // DEBUG
            testingClassifier(variationPerPos, pos);
            // DEBUG
            // TODO: Do anonymization on current values
            variationPerPos.remove(pos-SLIDING_WINDOW_LIMIT);
            if (p==SLIDING_WINDOW_LIMIT) p = 0;
            p++;
        }
    }

    @Override
    public Map<String, AnonymizedReadContainer> getAnonymizedReadContainer() {
        return anonReadContainer;
    }

    // DEBUG
    public void testingClassifier(Map<Integer, List<CalledVariation>> variationPerPos, int pos){
        List<CalledVariation> variationInPos = variationPerPos.get(pos);
        // System.out.println("pos " + pos);
        // if (variationInPos == null) {
        //     System.out.println("#POS_BUG: " + pos);
        // }
        for (CalledVariation var : variationInPos){
            System.out.print(var.toString());
            System.out.print("\tpresent in reads:\t");
            for (Map.Entry<String, Integer> entry : var.getSupportingReads().entrySet()){
                System.out.println(entry.getKey() + " in_read_pos=" + entry.getValue());
            }
        }
    }
    // DEBUG
}
