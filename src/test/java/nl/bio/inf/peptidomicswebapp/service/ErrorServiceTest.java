package nl.bio.inf.peptidomicswebapp.service;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;

import static org.junit.jupiter.api.Assertions.*;

/**
 * ErrorService Test class
 * @author Wouter Zeevat
 */
@SpringBootTest
class ErrorServiceTest {

    @Autowired
    private ErrorService errorService;

    @ParameterizedTest
    @CsvSource({"302,There was a code 302 error: \"Found\"", "304,There was a code 304 error: Moved Permanently", "404,There was a code 404 error: \"Not Found\". The resource you requested could not be found. Did you type its location correctly?"})
    void messagetest(String parms, String expected){
        String result = errorService.generateErrorMessage(Integer.parseInt(parms));
        assertEquals(result, expected);
    }
}