package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import nl.bio.inf.peptidomicswebapp.models.Plot;
import nl.bio.inf.peptidomicswebapp.service.PythonRunner;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.MediaType;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.RestController;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

@RestController
public class ResultController {


    @PostMapping(value = "/create_temp_file" , produces = MediaType.APPLICATION_JSON_VALUE)
    public void createTempFile(HttpServletRequest request, HttpSession session) throws IOException, InterruptedException {
        PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
        Path tempFilePath = Files.createTempFile(null, ".pdb");
        String tempUniqueName = String.valueOf(tempFilePath);
        tempUniqueName = tempUniqueName.substring(0, tempUniqueName.lastIndexOf(".pdb"));
        FileOutputStream fos = new FileOutputStream(tempFilePath.toFile());
        fos.write(pdb.getBytes());
        fos.close();
        File folderScripts = new ClassPathResource("scripts").getFile();
        File fullPath = null;
        for (File f: folderScripts.listFiles()) {
            if("read_pdb.py".equals(f.getName())) {
                fullPath = f;
            }
        }
        PythonRunner pythonRunner = new PythonRunner(fullPath.toString(), tempUniqueName, tempFilePath.toString());
        pythonRunner.startJobWithoutOutPut();
        Files.delete(tempFilePath);
        session.setAttribute("temp_numpyFile", tempUniqueName + ".npy");
    }

    @PostMapping(value = "/create_pca_plot")
    public @ResponseBody Plot createPca(HttpServletRequest request) throws IOException, InterruptedException {
        File folderScripts = new ClassPathResource("scripts").getFile();
        File fullPath = null;
        for (File f: folderScripts.listFiles()) {
            if("pca_dim_plot.py".equals(f.getName())) {
                fullPath = f;
            }
        }
        String numpyPath = request.getSession().getAttribute("temp_numpyFile").toString();
        PythonRunner pythonRunner = new PythonRunner(fullPath.toString(), numpyPath, "");
        String bytes = pythonRunner.startJobWithOutPut();
        return new Plot(bytes);
    }
}
