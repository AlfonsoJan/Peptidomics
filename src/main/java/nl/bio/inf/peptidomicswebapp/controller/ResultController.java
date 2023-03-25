package nl.bio.inf.peptidomicswebapp.controller;

import jakarta.servlet.http.HttpServletRequest;
import jakarta.servlet.http.HttpSession;
import nl.bio.inf.peptidomicswebapp.PeptidomicsWebAppApplication;
import nl.bio.inf.peptidomicswebapp.models.Chain;
import nl.bio.inf.peptidomicswebapp.models.PDB;
import nl.bio.inf.peptidomicswebapp.models.Plot;
import nl.bio.inf.peptidomicswebapp.service.PDBParser;
import nl.bio.inf.peptidomicswebapp.service.PythonService;
import org.springframework.beans.factory.annotation.Autowired;
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
import java.util.logging.Logger;

@RestController
public class ResultController {
    private static final Logger LOGGER  = Logger.getLogger(PeptidomicsWebAppApplication.class.getName());

    @Autowired
    private PythonService pythonService;

    @Autowired
    private PDBParser pdbParser;

    @PostMapping(value = "/get_stats_pdb")
    public PDB retrieveStatsPDB(HttpServletRequest request) {
        PDB pdb;
        try {
            pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            pdbParser.setParams(pdb);
            pdbParser.startFile();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return pdb;
    }

    @PostMapping(value = "/get_chains")
    public @ResponseBody Plot getChains(HttpServletRequest request){
        try {
            PDB pdb = (PDB) request.getSession().getAttribute("PDBFiles");
            File folderScripts  = new ClassPathResource("scripts").getFile();
            File fullPath = null;
            for (File f: folderScripts.listFiles()) {
                if("retrieve_chains_pdb.py".equals(f.getName())) {
                    fullPath = f;
                }
            }
            String chain = pythonService.getChainsPBD(String.valueOf(fullPath), pdb.getStructureId());
            return new Plot(chain);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }


    @PostMapping(value = "/create_temp_file" , produces = MediaType.APPLICATION_JSON_VALUE)
    public void createTempFile(HttpServletRequest request, HttpSession session) {
        try {
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
            pythonService.createTempNumpyFile(
                    fullPath.toString(),
                    tempUniqueName,
                    tempFilePath.toString(),
                    request.getSession().getAttribute("parameter").toString());
            Files.delete(tempFilePath);
            session.setAttribute("temp_numpyFile", tempUniqueName + ".npy");
        } catch (IOException ex) {
            LOGGER.severe("Error while creating a temp file, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    @PostMapping(value = "/create_pca_plot")
    public @ResponseBody Plot createPca(HttpServletRequest request) {
        try {
            File folderScripts = new ClassPathResource("scripts").getFile();
            File fullPath = null;
            for (File f: folderScripts.listFiles()) {
                if("pca_dim_plot.py".equals(f.getName())) {
                    fullPath = f;
                }
            }
            String numpyPath = request.getSession().getAttribute("temp_numpyFile").toString();
            String bytes = pythonService.createPcaPlot(fullPath.toString(), numpyPath);
            return new Plot(bytes);
        } catch (IOException ex) {
            LOGGER.warning("Error while reading creating pca plot, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    @PostMapping(value = "/create_scatter_plot")
    public @ResponseBody Plot createScatter(HttpServletRequest request) {
        try {
            File folderScripts = new ClassPathResource("scripts").getFile();
            File fullPath = null;
            for (File f: folderScripts.listFiles()) {
                if("axis_scatter_plot.py".equals(f.getName())) {
                    fullPath = f;
                }
            }
            String numpyPath = request.getSession().getAttribute("temp_numpyFile").toString();
            String bytes = pythonService.createScatterPlot(fullPath.toString(), numpyPath);
            return new Plot(bytes);
        } catch (IOException ex) {
            LOGGER.warning("Error while reading creating scatter plot, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    @PostMapping(value = "/pca_plotly_plot")
    public @ResponseBody Plot createPlotlyPcaPlot(HttpServletRequest request) {
        try {
            File folderScripts = new ClassPathResource("scripts").getFile();
            File fullPath = null;
            for (File f: folderScripts.listFiles()) {
                if("plotly_pca_plot.py".equals(f.getName())) {
                    fullPath = f;
                }
            }
            String numpyPath = request.getSession().getAttribute("temp_numpyFile").toString();
            String bytes = pythonService.createPlotlyPcaPlot(fullPath.toString(), numpyPath);
            return new Plot(bytes);
        } catch (IOException ex) {
            LOGGER.warning("Error while reading creating PCA plot, message=" + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }
}
