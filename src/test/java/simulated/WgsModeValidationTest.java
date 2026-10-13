package simulated;

/**
 * Simulated validation of the default mode of GenomeAnonymizer (sample type WGS)
 * @author Nicolas Gaitan
 */
public class WgsModeValidationTest extends AbstractWellformedValidation {

    @Override
    protected String getMode() {
        return "wgs";
    }

    @Override
    protected String[] getExtraArgs() {
        return new String[0];
    }
}
