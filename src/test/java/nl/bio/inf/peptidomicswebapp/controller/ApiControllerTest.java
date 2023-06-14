package nl.bio.inf.peptidomicswebapp.controller;

import nl.bio.inf.peptidomicswebapp.exceptions.EigenVectorsNotFoundException;
import nl.bio.inf.peptidomicswebapp.models.EigenVectors;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * ApiController Test class
 * @author Wouter Zeevat
 */
@SpringBootTest
public class ApiControllerTest {

    @Autowired
    ApiController controller;

    @Test
    @DisplayName("Tests a normal usage of the api controller function")
    void testApiController() {
        EigenVectors vector = controller.getEigenVector(1);
        assert(vector.getLength() == 1);
    }

    @ParameterizedTest
    @ValueSource(ints = {-1, 31, 60, 0, -100})
    @DisplayName("Tests the get file name function!")
    void testApiControllerInvalid(int integer) {
        assertThrows(EigenVectorsNotFoundException.class, () -> controller.getEigenVector(integer));
    }

    @ParameterizedTest
    @CsvSource({"0.1663496513136032,1", "0.0049478314358955,2", "-0.0003536110912536,3"})
    void testApiNormal(String expected, String params) {
        assert(Double.parseDouble(expected)== Double.parseDouble(controller.getEigenVector(Integer.valueOf(params)).getX().get(0)));
    }
}
