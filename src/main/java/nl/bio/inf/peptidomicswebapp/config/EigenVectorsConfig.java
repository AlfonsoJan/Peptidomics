package nl.bio.inf.peptidomicswebapp.config;

import nl.bio.inf.peptidomicswebapp.models.EigenVectors;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

/**
 *  This class read the eigen vectors once on startup and set it in a map
 * @author Jan Alfonso
 */

@Configuration
public class EigenVectorsConfig {

    @Value("${python.executable.folder}")
    private String pythonFolder;

    /**
     * This bean reads the files and set the eigenvectors in a map
     * @return Map with eigen vectors
     */
    @Bean
    public Map<Integer, EigenVectors> EigenVectorsMap() {
        Map<Integer, EigenVectors> eigenVectorsMap = new HashMap<>();
        Path filePath = Paths.get(pythonFolder, "vectors");
        try {
            for (int i = 1; i <= 30; i++) {
                Path vectorPath = Paths.get(filePath.toString(), "vectors_" + i + ".csv");
                EigenVectors eigenVectors = new EigenVectors(i);
                BufferedReader br = new BufferedReader(new FileReader(vectorPath.toString()));
                String line;
                while((line = br.readLine()) != null){
                    if (!(line.startsWith("#") | line.startsWith("X,Y,Z"))) {
                        eigenVectors.addLine(line);
                    }
                }
                eigenVectorsMap.put(i, eigenVectors);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return eigenVectorsMap;
    }
}
