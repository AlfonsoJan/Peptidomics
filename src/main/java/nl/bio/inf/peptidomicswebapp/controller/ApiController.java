package nl.bio.inf.peptidomicswebapp.controller;

import nl.bio.inf.peptidomicswebapp.exceptions.EigenVectorsNotFoundException;
import nl.bio.inf.peptidomicswebapp.models.EigenVectors;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RestController;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 *  This rest controller controls the api for the eigenvectors
 * @author Jan Alfonso Busker
 */

@RestController
public class ApiController {

    private final Map<Integer, EigenVectors> EigenVectorsMap;

    public ApiController(Map<Integer, EigenVectors> EigenVectorsMap) {
        this.EigenVectorsMap = EigenVectorsMap;
    }

    /**
     * Return all the eigenvectors
     * @return List of all the eigenvectors
     */
    @GetMapping("/api/v1/eigenvectors")
    List<EigenVectors> allEigenVectors() {
        return new ArrayList<>(EigenVectorsMap.values());
    }

    /**
     * Get the eigenvectors of a specific length or throws an error if its not found
     * @return EigenVector
     */
    @GetMapping("/api/v1/eigenvectors/{length}")
    EigenVectors getEigenVector(@PathVariable Integer length) {
        if ( EigenVectorsMap.containsKey(length) ) {
            return EigenVectorsMap.get(length);
        } throw new EigenVectorsNotFoundException(length);
    }
}
