package nl.bio.inf.peptidomicswebapp.service;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.stream.Stream;

public class PythonRunner implements CommandConstructor {
    private final String program = "python3";
    private final String options = "-u";
    private final String pdbPath;

    private final String pythonPath;
    private final String numpyPath;

    private String parameter;

    public PythonRunner(String pythonPath, String numpyPath, String pdbPath, String parameter) {
        this.pythonPath = pythonPath;
        this.numpyPath = numpyPath;
        this.pdbPath = pdbPath;
        this.parameter = parameter;
    }
    public PythonRunner(String pythonPath, String numpyPath, String pdbPath) {
        this.pythonPath = pythonPath;
        this.numpyPath = numpyPath;
        this.pdbPath = pdbPath;
    }

    private String[] constructCommandPrefix(){
        if (numpyPath.equals("")) {
            return new String[]{program, options, pythonPath};
        }
        return new String[]{program, options, pythonPath, numpyPath};
    }

    @Override
    public String[] constructCommand(){
        System.out.println(Arrays.toString(constructCommandPrefix()));
        return null;
    }

    public String startJobWithOutPut() throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder()
                .command(program, options, pythonPath, numpyPath);
        pb.redirectError(new File(numpyPath.substring(0, numpyPath.lastIndexOf('.')) + ".log"));
        Process p = pb.start();
        BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
        StringBuilder buffer = new StringBuilder();
        String line;
        while ((line = in.readLine()) != null){
            buffer.append(line);
        }
        int exitCode = p.waitFor();
        in.close();
        return buffer.toString();
    }

    public void startJobWithoutOutPut() throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder()
                .command(program, options, pythonPath, pdbPath, numpyPath, parameter);
        Process p = pb.start();
        int exitCode = p.waitFor();
    }
}
