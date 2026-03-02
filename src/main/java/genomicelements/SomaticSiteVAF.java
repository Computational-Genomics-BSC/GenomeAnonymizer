package genomicelements;

// Stores original VAF data and ALT base quality for a potential somatic site, used by fixVAF correction
public record SomaticSiteVAF(byte altAllele, float normalVAF, float tumorVAF,
                              byte medianNormalAltQuality, byte medianTumorAltQuality) {}
