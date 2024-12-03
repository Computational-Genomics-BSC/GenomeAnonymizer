package io;

import genomicelements.GenomicRegion;
import genomicelements.GenomicRegionBaseImpl;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Class that reads a bed file with genomic regions into GenomicRegion objects
 * @author Nicolas Gaitan
 */
public class GenomicRegionBedReader {

    public final static String BED_FIELD_SPARATOR = "\t";

    public static List<GenomicRegion> readGenomicRegionBED(String fileName) throws IOException {
        List<GenomicRegion> genomicRegions = new ArrayList<>();
        try(FileReader fileReader = new FileReader(fileName);
            BufferedReader reader = new BufferedReader(fileReader)){
            String line = reader.readLine();
            while(line != null){
                String[] elements = line.split(BED_FIELD_SPARATOR);
                String sequenceName = elements[0];
                int start =  Integer.parseInt(elements[1]);
                int end =  Integer.parseInt(elements[2]);
                GenomicRegion currentRegion = new GenomicRegionBaseImpl(sequenceName, start, end);
                genomicRegions.add(currentRegion);
                line = reader.readLine();
            }
        }
        return genomicRegions;
    }
}
