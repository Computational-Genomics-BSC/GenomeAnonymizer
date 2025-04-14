package genomicelements;

import htsjdk.samtools.Cigar;

public record PairInfoToUpdate(String readName, int alnStart, Cigar newCigar) { }
