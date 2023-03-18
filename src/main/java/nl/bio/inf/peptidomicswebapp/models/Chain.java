package nl.bio.inf.peptidomicswebapp.models;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

public class Chain {
    private final String chainId;
    private List<List<BigDecimal>> atoms = new ArrayList<>();
    private String[] seqres = new String[]{};
    private Long atomCount = 0L;


    public Chain(String chainId) {
        this.chainId = chainId;
    }

    public List<List<BigDecimal>> getAtoms() {
        return new ArrayList<>(atoms);
    }

    public void setAtoms(List<BigDecimal> atom) {
        this.atoms.add(atom);
    }

    public String getChainId() {
        return chainId;
    }

    public String[] getSeqres() {
        return seqres;
    }

    public void setSeqres(String[] seqresList) {
        this.seqres = concat(seqres, seqresList);
    }

    static String[] concat(String[]... arrays) {
        int length = 0;
        for (String[] array : arrays) {
            length += array.length;
        }
        String[] result = new String[length];
        int pos = 0;
        for (String[] array : arrays) {
            for (String element : array) {
                result[pos] = element;
                pos++;
            }
        }
        return result;
    }

    public void setCount() {
        this.atomCount = this.atomCount + 1;
    }

    public Long getCount() {
        return atomCount;
    }

    @Override
    public String toString() {
        return "Chain{" +
                "chainId='" +
                chainId +
                "', SEQRES Length='" +
                seqres.length +
                "', ATOM length=" +
                atomCount +
                "residues'}";
    }
}
