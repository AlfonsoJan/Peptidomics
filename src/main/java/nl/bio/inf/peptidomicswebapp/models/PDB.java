package nl.bio.inf.peptidomicswebapp.models;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;

public class PDB {
    private static final String DOWNLOAD_BY_ID_URL = "https://files.rcsb.org/download/%s.pdb";

    private String structureId;

    public PDB(String structureId) {
        this.structureId = structureId;
    }

    public InputStream getInputStream() throws IOException {
        URL url = new URL(String.format(DOWNLOAD_BY_ID_URL, this.structureId));
        URLConnection connection = url.openConnection();
        return connection.getInputStream();
    }

    public byte[] getBytes() throws IOException {
        return getInputStream().readAllBytes();
    }

    public String getStructureId() {
        return this.structureId;
    }
}
