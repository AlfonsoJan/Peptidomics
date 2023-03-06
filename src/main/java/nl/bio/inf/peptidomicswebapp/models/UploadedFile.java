package nl.bio.inf.peptidomicswebapp.models;

public record UploadedFile(String originalFilename, byte[] data) {

    public byte[] getBytes() {
        return data.clone();
    }

}

//public class UploadedFile {
//    private final String originalFilename;
//    private final byte[] data;
//
//    public UploadedFile(String fileName, byte[] data) {
//        this.originalFilename = fileName;
//        this.data = data.clone();
//    }
//
//    public String getOriginalFilename() {
//        return originalFilename;
//    }
//
//    public byte[] getData() {
//        return data.clone();
//    }
//
//    public int getSize() {
//        return data.length;
//    }
//}