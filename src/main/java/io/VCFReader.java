package io;

import genomicelements.CalledVariation;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class VCFReader {
    public static Map<String, Map<Integer, CalledVariation>> readVCF(String vcfFilePath)throws IOException {
        Map<String, Map<Integer, CalledVariation>> somaticCallsToKeep = new HashMap<>();
        File vcfFile = new File(vcfFilePath);
        try(VCFFileReader reader = new VCFFileReader(vcfFile, false)){
            for(VariantContext variantContext : reader){
                CalledVariation calledVar = new CalledVariation(variantContext);
                Map<Integer, CalledVariation> perPosMap = somaticCallsToKeep.computeIfAbsent(calledVar.getSeqName(), v -> new HashMap<>());
                perPosMap.put(calledVar.getPos(), calledVar);
                //DEBUG
                //System.out.println("& Variant=" + calledVar);
                //DEBUG
            }
        }
        return somaticCallsToKeep;
    }
}
