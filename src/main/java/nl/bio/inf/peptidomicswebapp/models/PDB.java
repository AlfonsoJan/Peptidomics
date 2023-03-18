package nl.bio.inf.peptidomicswebapp.models;

import io.micrometer.common.util.StringUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class PDB {
    private static final String DOWNLOAD_BY_ID_URL = "https://files.rcsb.org/download/%s.pdb";

    private final String structureId;
    private final byte[] bytes;
    private final String fileName;
    private List<Chain> chainList = new ArrayList<>();

    public PDB(String structureId) throws IOException {
        this.structureId = structureId;
        this.bytes = getBytesConnection();
        this.fileName = structureId + ".pdb";
    }

    public PDB(String fileName, byte[] fileBytes, String structureId){
        this.fileName = fileName;
        this.bytes  = fileBytes;
        this.structureId = structureId;
    }

    public void setChainList(List<Chain> chainList) {
        this.chainList = chainList;
    }

    public List<Chain> getChainList() {
        return chainList;
    }

    public static String getStructureFromInputstream(InputStream is) throws IOException {
        String line;
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(is));
        while( (line = bufferedReader.readLine()) != null ) {
            if (line.startsWith("HEADER")) {
                return line.strip().substring(line.strip().lastIndexOf(" ")+1);
            }
        }
        return "";
    }

    private InputStream getInputStream() throws IOException {
        URL url = new URL(String.format(DOWNLOAD_BY_ID_URL, this.structureId));
        URLConnection connection = url.openConnection();
        return connection.getInputStream();
    }

    private byte[] getBytesConnection() throws IOException {
        return getInputStream().readAllBytes();
    }

    public byte[] getBytes(){
        return bytes;
    }

    public String getFileName() {
        return fileName;
    }

    public String getStructureId() {
        return this.structureId;
    }
}
