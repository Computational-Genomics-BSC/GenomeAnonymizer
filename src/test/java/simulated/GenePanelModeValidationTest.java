package simulated;

/**
 * Simulated validation of the gene panel mode of GenomeAnonymizer with VAF correction, using the command that the
 * README recommends for gene panels. The dataset has WGS depth, so this does not evaluate how well the mode works
 * on real panels: it records its behaviour to detect changes
 * @author Nicolas Gaitan
 */
public class GenePanelModeValidationTest extends AbstractWellformedValidation {

    @Override
    protected String getMode() {
        return "gene_panel";
    }

    @Override
    protected String[] getExtraArgs() {
        return new String[]{"-sampleType", "gene_panel", "-maxDepth", "-1", "-fixVAF"};
    }
}
