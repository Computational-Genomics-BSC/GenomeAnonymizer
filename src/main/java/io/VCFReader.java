package io;

import genomicelements.CalledVariation;
import genomicelements.CalledVariation.VariantType;
import genomicelements.CalledVariation.BreakendSVRecord;
import genomicelements.CalledVariation.ShorthandSVRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.IllegalFormatException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * VCFReader that comprehensively parses Standard, Breakend and shorthand VCF records into CalledVariation objects
 * @author Nicolas Gaitan
 * @author Rodrigo Martin
 */
public class VCFReader {

    // Regex for SVs
    private static final Pattern BREAKEND_SV_REGEX = Pattern.compile("([.A-Za-z]*)(\\[|\\])([^\\]\\[:]+:[0-9]+)(\\[|\\])([.A-Za-z]*)");
    private static final Pattern SHORTHAND_SV_REGEX = Pattern.compile("<(DEL|INS|DUP|INV|CNV)(:[A-Za-z0-9]+)*>");
    private static final Pattern SGL_SV_REGEX = Pattern.compile("\\.[.A-Za-z]+|[.A-Za-z]+\\.");
    private static final Pattern STANDARD_RECORD_REGEX = Pattern.compile("([.A-Za-z]+)");

    public static Map<String, Map<Integer, CalledVariation>> readVCF(String vcfFilePath)
            throws Exception {
        Map<String, Map<Integer, CalledVariation>> somaticCallsToKeep = new HashMap<>();
        File vcfFile = new File(vcfFilePath);
        try(VCFFileReader reader = new VCFFileReader(vcfFile, false)){
            int i = 0;
            for(VariantContext variantContext : reader){
                CalledVariation calledVar;
                if (variantContext.isSNP() || variantContext.isIndel()){
                    calledVar = processStandardVariant(variantContext);
                }
                else{
                    calledVar = processBreakendSV(variantContext);
                    if (calledVar == null) calledVar = processShorthandSV(variantContext);
                    if (calledVar == null) calledVar = processSGLSV(variantContext);
                    if (calledVar == null) throw new ParseException("Complex Variant record in vcf is not parseable: "
                            + variantContext.getContig() + " " + variantContext.getAlternateAllele(0).getDisplayString() + "\n"
                            + "Displaying line number fo malformed record: ", i);
                }
                // By default, redundant records are saved as one CalledVariation only (Last one is saved)
                Map<Integer, CalledVariation> perPosMap = somaticCallsToKeep.computeIfAbsent(calledVar.getSeqName(), v -> new HashMap<>());
                perPosMap.put(calledVar.getPos(), calledVar);
                i++;
            }
        }
        return somaticCallsToKeep;
    }

    /**
     * Method to create a CalledVarition object from the fields in a Standard VariantContext (SNV|INDEL) from a VCF file
     * @param varContext
     */
    public static CalledVariation processStandardVariant(VariantContext varContext){
        CalledVariation answer = new CalledVariation(varContext.getContig(), varContext.getStart(), varContext.getEnd(),
                VariantType.SNV, varContext.getLengthOnReference(),
                varContext.getAlternateAllele(0).getBases(), new byte[0]);
        answer.setVariantType(getTypeFromVarContext(varContext));
        if(VariantType.INS == answer.getVariantType()) {
            answer.setEnd(answer.getEnd()+1);
            answer.setLength(answer.getAllele().length - 1);
        }
        answer.setSomaticVariationType(CalledVariation.SomaticVariationType.NOT_SOMATIC);
        answer.setDiffusedStatus(false);
        answer.setLinkToAnotherGermline(false);
        return answer;
    }

    private static VariantType getTypeFromVarContext(VariantContext context){
        VariantType varType = VariantType.SNV;
        if(context.isSNP()) varType = VariantType.SNV;
        if(context.isSimpleDeletion()) varType = VariantType.DEL;
        if(context.isSimpleInsertion()) varType = VariantType.INS;
        if(context.isSymbolicOrSV()) varType = VariantType.SV;
        if(context.isComplexIndel()){
            int length;
            int altLength = context.getAlternateAllele(0).length();
            int refLength = context.getReference().length();
            if(altLength > refLength){
                varType = VariantType.INS;
            }
            else{
                varType = VariantType.DEL;
            }
        }
        return varType;
    }

    /**
     * Method to create a CalledVarition object from parsing a Breakend VariantContext record (SV) from a VCF file
     * @param varContext
     */
    public static CalledVariation processBreakendSV(VariantContext varContext) {
        Matcher svMatchBreakend = BREAKEND_SV_REGEX.matcher(varContext.getAlternateAllele(0).getDisplayString());
        if (!svMatchBreakend.matches()) {
            return null;
        }
        // Extract ALT data from regex
        String altPrefix = svMatchBreakend.group(1);
        String altBracket = svMatchBreakend.group(2);
        String [] altContigAndPos = svMatchBreakend.group(3).split(":");
        String altContig = altContigAndPos[0];
        int altPos = Integer.parseInt(altContigAndPos[1]);
        String altSuffix = svMatchBreakend.group(5);
        BreakendSVRecord altSvBreakend = new BreakendSVRecord(altPrefix, altBracket, altContig,
                altPos, altSuffix);
        // End position
        int endPos = altPos;
        // Extract type
        VariantType variantType;
        int length;
        if (!varContext.getContig().equals(altContig)) {
            // BND & INS with different contig ~ TRA
            variantType = VariantType.TRA;
            length = 0;
        } else {
            // INV -> 1 10 N]1:20] or 1 20 N]1:10]
            //        1 10 [1:20[N or 1 20 [1:10[N
            // DEL -> 1 10 N[1:20[ or 1 20 ]1:10]N
            // DUP -> 1 10 ]1:20]N or 1 20 N[1:10[
            length = Math.abs(endPos - varContext.getStart());
            String equivalentBracket = altBracket;
            String equivalentPrefix = altPrefix;
            // Transform REF/ALT to equivalent notation
            if (altPos < varContext.getStart()) {
                equivalentBracket = altBracket.equals("[") ? "]" : "[";
                equivalentPrefix = altSuffix;
            }
            if (!equivalentPrefix.isEmpty() && equivalentBracket.equals("[")) {
                variantType = VariantType.DEL;
            } else if (equivalentPrefix.isEmpty() && equivalentBracket.equals("]")) {
                variantType = VariantType.DUP;
            } else {
                variantType = VariantType.INV;
            }
        }
        // Create new record
        CalledVariation answer = new CalledVariation(varContext.getContig(), varContext.getStart(), endPos, variantType, length,
                varContext.getAlternateAllele(0).getDisplayBases(), varContext.getReference().getDisplayBases());
        answer.setBreakendRecord(altSvBreakend);
        return answer;
    }

    /**
     * Method to create a CalledVarition object from parsing a Shorthand notation VariantContext record (SV) from a VCF file
     * @param varContext
     */
    public static CalledVariation processShorthandSV(VariantContext varContext) {
        Matcher svMatchShorthand = SHORTHAND_SV_REGEX.matcher(varContext.getAlternateAllele(0).getDisplayString());
        if (!svMatchShorthand.matches()) {
            return null;
        }
        // Extract ALT data from regex
        String altType = svMatchShorthand.group(1);
        String[] altExtra = svMatchShorthand.group(2) != null ?
                svMatchShorthand.group(2).substring(1).split(":") : null;
        ShorthandSVRecord altSvShorthand = new ShorthandSVRecord(altType, altExtra);
        // TODO: Check with a VCF that has this notation if VariantContext gives the propper corodinates
        int length = Math.abs(varContext.getEnd() - varContext.getStart());
        // Extract type
        VariantType variantType = switch (altType) {
            case "DEL" -> VariantType.DEL;
            case "INS" -> {
                if (!varContext.hasAttribute("SVLEN")) {
                    System.err.println("Warning: SVLEN not found in INFO field for <INS> shorthand record. Defaults to 0.");
                    length = 0;
                } else {
                    int svlen = varContext.getAttributeAsInt("SVLEN", 0);
                    length = Math.abs(svlen);
                }
                yield VariantType.INS;
            }
            case "DUP" -> VariantType.DUP;
            case "INV" -> VariantType.INV;
            case "CNV" -> VariantType.CNV;
            default -> throw new IllegalArgumentException("Unknown variant type: " + altType + ". Skipping:\n" + varContext.toString());
        };
        // Create new record
        CalledVariation answer = new CalledVariation(varContext.getContig(), varContext.getStart(),
                varContext.getEnd(), variantType, length, varContext.getAlternateAllele(0).getDisplayBases(),
                varContext.getReference().getDisplayBases());
        answer.setShortHandRecord(altSvShorthand);
        return answer;
    }

    /**
     * Method to create a CalledVarition object from parsing a Single Breakpoint VariantContext record (SV) from a VCF file
     * @param varContext
     */
    public static CalledVariation processSGLSV(VariantContext varContext) {
        Matcher svMatchSgl = SGL_SV_REGEX.matcher(varContext.getAlternateAllele(0).getDisplayString());
        if (!svMatchSgl.matches() || !varContext.hasAttribute("SVTYPE")) {
            return null;
        }
        VariantType variantType = VariantType.SGL;
        int length = 0;
        // Create new record
        return new CalledVariation(varContext.getContig(), varContext.getStart(), varContext.getEnd(),
                variantType, length, varContext.getAlternateAllele(0).getDisplayBases(),
                varContext.getReference().getDisplayBases());
    }
}
