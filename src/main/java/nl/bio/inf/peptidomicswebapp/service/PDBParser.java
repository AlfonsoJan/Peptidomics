package nl.bio.inf.peptidomicswebapp.service;

import nl.bio.inf.peptidomicswebapp.models.Chain;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import org.springframework.stereotype.Service;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;

@Service
public class PDBParser {
    private PDB pdb;
    private Map<String, Integer> map;
    private ArrayList<Chain> chains;

    public void setParams(PDB pdb) {
        this.pdb = pdb;
        this.map = new HashMap<>();
        this.chains = new ArrayList<>();
    }

    public void startFile() throws IOException {
        InputStream input = new ByteArrayInputStream(pdb.getBytes());
        BufferedReader reader = new BufferedReader(new InputStreamReader(input));
        String line;
        while ((line = reader.readLine()) != null) {
            parseLine(line);
        }
        reader.close();
        pdb.setChainList(chains);
    }

    private void parseLine(String line) {
        if (line.startsWith("COMPND")) {
            parseHeader(line);
        } else if (line.toUpperCase().startsWith("SEQRES")) {
            parseSeqres(line);
        } else if (line.toUpperCase().startsWith("ATOM")) {
            parseAtom(line);
        }
    }

    private void parseAtom(String line) {
        String[] elem = line.split("\\s+");
        if (elem.length == 10) {
            parseAtomlength10(line);
        } else if (elem.length == 11) {
            parseAtomlength11(line);
        } else if (elem.length == 12) {
            parseAtomlength12(line);
        }
    }

    private void parseAtomlength12(String line) {
        String[] elem = line.split("\\s+");
        String chain = elem[4];
        int current = map.get(chain);
        boolean ifDifferent = false;
        if (current != Integer.parseInt(elem[5])) {
            ifDifferent = true;
            map.put(chain, Integer.parseInt(elem[5]));
        }
        for (Chain c: chains ) {
            if (c.getChainId().equals(chain)) {
                if (ifDifferent) {
                    c.setCount();
                }
                c.setAtoms(List.of(
                        new BigDecimal(elem[6]),
                        new BigDecimal(elem[7]),
                        new BigDecimal(elem[8])
                ));
            }
        }
    }

    private void parseAtomlength11(String line) {
        String[] elem = line.split("\\s+");
        String chain = elem[4];
        if (elem[4].strip().length() > 1) {
            chain = elem[4].substring(0, 1);
        }
        int current = map.get(chain);
        boolean ifDifferent = getSequenceNumberAtom11(elem, current, chain);
        getCoordinatesAtom11(elem, ifDifferent, chain);
    }

    private void getCoordinatesAtom11(String[] elem, boolean ifDifferent, String chain) {
        try {
            Float.parseFloat(elem[elem.length - 2]);
            for (Chain c: chains ) {
                if (c.getChainId().equals(chain)) {
                    if (ifDifferent) {
                        c.setCount();
                    }
                    c.setAtoms(List.of(
                            new BigDecimal(elem[5]),
                            new BigDecimal(elem[6]),
                            new BigDecimal(elem[7])
                    ));
                }
            }
        } catch (Exception e) {
            for (Chain c: chains ) {
                if (c.getChainId().equals(chain)) {
                    if (ifDifferent) {
                        c.setCount();
                    }
                    c.setAtoms(List.of(
                            new BigDecimal(elem[6]),
                            new BigDecimal(elem[7]),
                            new BigDecimal(elem[8])
                    ));
                }
            }
        }
    }

    private boolean getSequenceNumberAtom11(String[] elem, int current, String chain) {
        try {
            if (current != Integer.parseInt(elem[5])) {
                map.put(chain, Integer.parseInt(elem[5]));
                return true;
            }
            return false;
        } catch (Exception e) {
            if (current != Integer.parseInt(elem[4].substring(1))) {
                map.put(chain, Integer.parseInt(elem[4].substring(1)));
                return true;
            }
            return false;
        }
    }

    private void parseAtomlength10(String line) {
        String[] elem = line.split("\\s+");
        String chain = elem[3];
        if (elem[3].strip().length() != 1) {
            chain = elem[4].substring(0, 1);
        }
        int current = map.get(chain);
        boolean ifDifferent = getSequenceNumberAtom10(elem, current, chain);
        getCoordinatesAtom10(elem, ifDifferent, chain);
    }

    private void getCoordinatesAtom10(String[] elem, boolean ifDifferent, String chain) {
        try {
            Float.parseFloat(elem[elem.length - 2]);
            for (Chain c: chains ) {
                if (c.getChainId().equals(chain)) {
                    if (ifDifferent) {
                        c.setCount();
                    }
                    c.setAtoms(List.of(
                            new BigDecimal(elem[4]),
                            new BigDecimal(elem[5]),
                            new BigDecimal(elem[6])
                    ));
                }
            }
        } catch (Exception e) {
            for (Chain c: chains ) {
                if (c.getChainId().equals(chain)) {
                    if (ifDifferent) {
                        c.setCount();
                    }
                    c.setAtoms(List.of(
                            new BigDecimal(elem[5]),
                            new BigDecimal(elem[6]),
                            new BigDecimal(elem[7])
                    ));
                }
            }
        }
    }

    private boolean getSequenceNumberAtom10(String[] elem, int current, String chain) {
        try {
            Float.parseFloat(elem[elem.length - 2]);
            if (current != Integer.parseInt(elem[4])) {
                map.put(chain, Integer.parseInt(elem[4]));
                return true;
            }
            return false;
        } catch (Exception e) {
            if (current != Integer.parseInt(elem[4].substring(1))) {
                map.put(chain, Integer.parseInt(elem[4].substring(1)));
                return true;
            }
            return false;
        }
    }

    private void parseSeqres(String line) {
        String[] elem = line.split("\\s+");
        String chain  = elem[2];
        for (Chain c: chains ) {
            if (c.getChainId().equals(chain)) {
                c.setSeqres(Arrays.copyOfRange(elem, 4, elem.length));
            }
        }
    }

    private void parseHeader(String line) {
        String[] elem = line.split("\\s+");
        boolean contains = Arrays.stream(elem).anyMatch("CHAIN:"::equals);
        if (contains) {
            line = line.substring(line.lastIndexOf(":") + 1).strip();
            for (String chain: line.split(",")) {
                chains.add(new Chain(chain.strip().replace(";", "")));
                map.put(chain.strip().replace(";", ""), 0);
            }
        }
    }
}