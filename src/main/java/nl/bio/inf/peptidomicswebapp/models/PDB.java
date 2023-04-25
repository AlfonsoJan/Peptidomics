package nl.bio.inf.peptidomicswebapp.models;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 *  This class handles the creation of the temporary pdb files.
 * @author Jan Alfonso Busker
 */
public class PDB {
    private static final String DOWNLOAD_BY_ID_URL = "https://files.rcsb.org/download/%s.pdb";

    private final String structureId;
    private final byte[] bytes;
    private final String fileName;

    /**
     * Constructor if a PDB code is given
     * @param structureId
     * @throws IOException
     */
    public PDB(String structureId) throws IOException {
        this.structureId = structureId;
        this.bytes = getBytesConnection();
        this.fileName = structureId + ".pdb";
    }

    /**
     * Constructor if a file is uploaded
     * @param fileBytes
     * @param fileName
     * @throws IOException
     */
    public PDB(byte[] fileBytes, String fileName) throws IOException {
        this.structureId = getStructureFromInputstream(fileBytes);
        this.bytes = fileBytes;
        this.fileName = fileName;
    }

    /**
     * Return the inputstream of the pdb file
     * @return
     * @throws IOException
     */
    private InputStream getInputStream() throws IOException {
        URL url = new URL(String.format(DOWNLOAD_BY_ID_URL, this.structureId));
        URLConnection connection = url.openConnection();
        return connection.getInputStream();
    }


    /**
     * Returns the bytes of the inputstream
     * @return
     * @throws IOException
     */
    private byte[] getBytesConnection() throws IOException {
        return getInputStream().readAllBytes();
    }


    /**
     * This method will create a temporary file, with a random prefix
     * @return
     */
    public String createTempFile() {
        try {
            Path tempFilePath = Files.createTempFile(null, ".pdb");
            FileOutputStream fos = new FileOutputStream(tempFilePath.toFile());
            BufferedReader reader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(bytes)));
            // Read the stream and write to the temporary file
            while(reader.ready()) {
                String line = reader.readLine();
                if (line.toUpperCase().startsWith("ATOM")) {
                    fos.write((line + "\n").getBytes(StandardCharsets.UTF_8));
                }
            }
            reader.close();
            fos.close();
            return tempFilePath.toString();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // Here below are all the static method needed for the constructor.

    /**
     * Returns the inputstream of the pdb file
     * @param pdbCode
     * @return
     * @throws IOException
     */
    private static byte[] getBytesConnection(String pdbCode) throws IOException {
        return getInputStream(pdbCode).readAllBytes();
    }

    /**
     * Returns the bytes of the inputstream
     * @param pdbCode
     * @return
     * @throws IOException
     */
    private static InputStream getInputStream(String pdbCode) throws IOException {
        URL url = new URL(String.format(DOWNLOAD_BY_ID_URL, pdbCode));
        URLConnection connection = url.openConnection();
        return connection.getInputStream();
    }

    /**
     * Get the structure id of the pdb file that is uploaded
     * @param fileBytes
     * @return
     * @throws IOException
     */
    public static String getStructureFromInputstream(byte[] fileBytes) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fileBytes)));
        String line;
        while( (line = reader.readLine()) != null ) {
            // The id is on the same line as the line that starts with header
            if (line.startsWith("HEADER")) {
                line = line.strip().substring(line.strip().lastIndexOf(" ")+1);
                break;
            }
        }
        reader.close();
        return line;
    }

    /**
     * The static method to create a file. Needed in createTemporaryFileCompare.
     * @param pdbCode
     * @return
     */
    public static String createTempFile(String pdbCode) {
        try {
            byte[] fileBytes = getBytesConnection(pdbCode);
            Path tempFilePath = Files.createTempFile(null, ".pdb");
            FileOutputStream fos = new FileOutputStream(tempFilePath.toFile());
            BufferedReader reader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fileBytes)));
            while(reader.ready()) {
                String line = reader.readLine();
                if (line.toUpperCase().startsWith("ATOM")) {
                    fos.write((line + "\n").getBytes(StandardCharsets.UTF_8));
                }
            }
            reader.close();
            fos.close();
            return tempFilePath.toString();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public String getStructureId() {
        return this.structureId;
    }

    public byte[] getBytes(){
        return bytes.clone();
    }

    public String getFileName() {
        return fileName;
    }
}