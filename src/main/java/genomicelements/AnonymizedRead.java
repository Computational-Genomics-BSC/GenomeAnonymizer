package genomicelements;

import java.util.List;
import java.util.Map;

// May be better as a parent class
public interface AnonymizedRead{

    public void setVariantsToAnonymize(List<CalledVariation> variants);

    //public void modifyBaseInRead(int inReadPosition, byte asciiBase);

    public void anonymizeVariants();

    public boolean isAnonymized();
}
