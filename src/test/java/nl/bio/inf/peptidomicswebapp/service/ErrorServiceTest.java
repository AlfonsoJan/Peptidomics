package nl.bio.inf.peptidomicswebapp.service;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.*;

/**
Wouter Zeevat
 **/
class ErrorServiceTest {

    @Disabled
    @ParameterizedTest
    @ValueSource(strings ={"302", "304", "404"})
    void generateErrorMessage(String parms) {
        String result = new ErrorService().generateErrorMessage(Integer.parseInt(parms));
        String message = "There was a code " + parms + " error: ";
        //assertEquals(result, message);
    }
}