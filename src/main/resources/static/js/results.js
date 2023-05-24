// author: Jan Alfonso Busker
//
// author: Wouter Zeevat
toastr.options = {
    "closeButton": true,
    "debug": false,
    "newestOnTop": true,
    "progressBar": false,
    "positionClass": "toast-bottom-right",
    "preventDuplicates": false,
    "onclick": null,
    "showDuration": "300",
    "hideDuration": "1000",
    "timeOut": "5000",
    "extendedTimeOut": "1000",
    "showEasing": "linear",
    "hideEasing": "linear",
    "showMethod": "fadeIn",
    "hideMethod": "fadeOut"
}
let resultJS = {
    aminoAcidCodes: [
        "ALA", // Alanine
        "ARG", // Arginine
        "ASN", // Asparagine
        "ASP", // Aspartic Acid
        "CYS", // Cysteine
        "GLN", // Glutamine
        "GLU", // Glutamic Acid
        "GLY", // Glycine
        "HIS", // Histidine
        "ILE", // Isoleucine
        "LEU", // Leucine
        "LYS", // Lysine
        "MET", // Methionine
        "PHE", // Phenylalanine
        "PRO", // Proline
        "SER", // Serine
        "THR", // Threonine
        "TRP", // Tryptophan
        "TYR", // Tyrosine
        "VAL", // Valine
    ],

    points: [],
    hidden: false,
    reset: true,
}

// Function that converts hue to RGB
function HSLToRGB (h, s, l) {
    s /= 100;
    l /= 100;
    const k = n => (n + h / 30) % 12;
    const a = s * Math.min(l, 1 - l);
    const f = n =>
        l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
    return [255 * f(0), 255 * f(8), 255 * f(4)];
};

// Function that creates perfectly distanced colors that follow the rainbow pattern
// according to an array of keys, creates {"A": red, "B": green, "C": pink}
function createColors(keys) {

    let colorArray = [];
    for (let i = 0; i < keys.length; i++) {
        let hsl_value = 255 / keys.length * i;
        let rgb = HSLToRGB(hsl_value, 100, 70);
        colorArray.push(rgb)
    }

    return keys.map((x, i) => ({x, y: colorArray[i]}));
}

// Selects a certain part of the JSMOL
function selectView(datapoints) {

    // Unhides the JMOL to rehide it after! BUG FIX
    let was_hidden = false;
    if (resultJS.hidden) {
        hideGlobal();
        was_hidden = true;
    }

    let a = document.createElement("a");
    let script = `"spacefill off; select all; wireframe 0.05; cartoons off; color [84,84,84]; `

    // Select each amino acid for all positions
    resultJS.points = datapoints;
    resultJS.points.forEach(point => {
        script += `select atomno>${parseInt(point[0]) - 1} and atomno<${parseInt(point[1]) + 1}; cartoons; color [10,0,255]; `;
    })
    script += "\"";
    a.href = `javascript:Jmol.script(jmol1, ${script})`
    a.click();
    document.getElementById("zoom-btn").classList.remove("is-hidden");
    document.getElementById("3d-text").innerText = `${resultJS.points.length} selected`;

    // Rehides
    if (was_hidden) {
        hideGlobal();
        script = `"zoom 0"`;
    }
    resultJS.reset = false;
}
// Function that create the 3d scatter plotly plot
function create3dPlot(result, colors, scores) {
    let elem = document.getElementById("spinner-scatter-3d");
    elem.parentNode.removeChild(elem);
    document.getElementById("placeholder-scatter-3d").style.display= '';
    let { data, initialView, secondaryView, thirdView } = getDataPCA3D(result, colors, scores);
    let updateMenus = initializePlotlyButtons(initialView, secondaryView, thirdView);
    let layout = setLayoutPlotly3D(updateMenus);
    let config = {responsive: true};
    Plotly.newPlot("placeholder-scatter-3d", data, layout, config);
    let myPlot = document.getElementById("placeholder-scatter-3d");

    // Click event
    myPlot.on("plotly_click", function(datapoints){
        if (datapoints === undefined) return;
        if (!datapoints.points[0].data.hasOwnProperty("freetext")) return; // If its a gray dot AKA background data
        selectView([datapoints.points[0].data.freetext[datapoints.points[0].pointNumber]]);

    });
}
// Function that create the 2d scatter plotly plot for the PCA results
function create2dPlot(result, colors, scores) {
    let elem = document.getElementById("spinner-scatter");
    elem.parentNode.removeChild(elem);
    document.getElementById("placeholder-scatter").style.display= '';
    let { data, initialView, secondaryView, thirdView } = getDataPCA2D(result, colors, scores);
    let updateMenus = initializePlotlyButtons(initialView, secondaryView, thirdView);
    let layout = setLayoutPlotly2D(updateMenus);
    let config = {responsive: true};
    Plotly.newPlot('placeholder-scatter', data, layout, config);
    let myPlot = document.getElementById("placeholder-scatter");

    // Click event
    myPlot.on("plotly_click", function(datapoints){
        console.log(datapoints)
        if (!datapoints.points[0].data.hasOwnProperty("freetext")) return; // If its a gray dot AKA background data
        selectView([datapoints.points[0].data.freetext[datapoints.points[0].pointNumber]])

    });

    // Select multiple with lasso or box select event
    myPlot.on("plotly_selected", function(datapoints){
        if (datapoints === undefined) return;
        if (datapoints.points.length < 1) {
            return;
        }

        let points = [];
        datapoints.points.forEach(data => {
            let index = data.pointIndex;
            if (!data.data.hasOwnProperty("freetext")) return; // If its a gray dot AKA background data
            let point = [data.data.freetext[index][0], data.data.freetext[index][1]]
            points.push(point);
        });
        selectView(points);
    });
}
// Function that set the chains on the site
function setChain(chains, colors) {
    // The size for how many chains per row
    const chunkSize = 4
    const chain = Object.keys(JSON.parse(chains["bytes"])).map((key) => [key, JSON.parse(chains["bytes"])[key]]);
    if (chain.length < 1) {
        let element = document.getElementById("stats-spinner");
        element.parentElement.removeChild(element);
        document.getElementById("place-text").style.display = '';
        return;
    }
    let elem = document.getElementById("pdb-stats");
    elem.parentElement.removeChild(elem);
    for (let i = 0; i < chain.length; i += chunkSize) {
        const chunk = chain.slice(i, i + chunkSize);
        const columns = document.createElement("div");
        columns.className = "columns";

        // Creates element for each chain
        chunk.forEach(function (c) {
            const column = document.createElement("div");
            column.className = "column";
            columns.appendChild(column);

            const card = document.createElement("div");
            card.className = "card clickable";
            card.onclick = clickChain;
            column.appendChild(card);

            const cardContent = document.createElement("div");
            cardContent.className = "chain-content card-content";
            card.appendChild(cardContent);

            // Loops rainbow colors to get best colors that the plot will late be using
            let rgb = colors.filter(word => word.x === c[0])[0]["y"];
            const chainId = document.createElement("p");
            chainId.className = "is-size-5 has-text-weight-bold";
            chainId.textContent = `Chain: ${c[0]}`;
            chainId.style.backgroundColor = `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
            cardContent.appendChild(chainId);

            const atomLength = document.createElement("p");
            atomLength.className = "is-size-5";
            atomLength.textContent = `ATOM: ${c[1]} residues`;
            cardContent.appendChild(atomLength);
        });
        document.getElementById("stats-pdb").appendChild(columns);
    }
}

function clickChain(item) {
    let ele = item.target;
    while (!ele.parentElement.classList.contains("card")) {
        ele = ele.parentElement;
    }
    let chain = ele.children[0].textContent.charAt(ele.children[0].textContent.length - 1);
    let color = ele.children[0].style.backgroundColor.replace("rgb(", "").replace(")", "").replace(" ", "");
    resetScript();

    let a = document.createElement("a");
    let script = `"select all; color [90,90,90]; select chain=${chain}; cartoons only; color [${color}]; background white; zoom 0"`;
    a.href = `javascript:Jmol.script(jmol1, ${script})`
    a.click();
    resultJS.reset = true;

}

// Resets the JSMOL view
function resetScript() {

    // Unhides if it was hidden to prevent bugs
    if (resultJS.hidden) {
        hideGlobal();
    }
    document.getElementById("zoom-btn").classList.add("is-hidden");
    let a = document.createElement("a");
    let script = `"select all; cartoons only; color structure; background white; zoom 0"`;
    a.href = `javascript:Jmol.script(jmol1, ${script})`
    a.click();
    document.getElementById("3d-text").innerText = "0 selected";
    resultJS.reset = true;
}

// Zooms in on the selected part and hides the other parts in JSMOL
function hideGlobal() {

    let a = document.createElement("a");

    // Don't go through if nothing is selected
    if (resultJS.reset) {
        return;
    }

    let script;
    if (!resultJS.hidden) {
        console.log("Hidden the global protein!")
        script = `"hide all; display `

        // Adds selected part each time because JSMOL does not allow you to do it once
        scriptParts = []
        resultJS.points.forEach(point => {
            scriptParts.push(`atomno>${parseInt(point[0]) - 1} and atomno<${parseInt(point[1]) + 1}`);
        });
        script += `${scriptParts.join(", ")}; color [10,0,255]"`;
        a.href = `javascript:Jmol.script(jmol1, ${script})`
        resultJS.hidden = true;
        a.click()

    } else {
        console.log("Unhidden the global protein!")
        resultJS.hidden = false;
        script = '"hide none"'
        a.href = `javascript:Jmol.script(jmol1, ${script})`
        a.click();
    }
    script = `"zoom 0"`;
    a.href = `javascript:Jmol.script(jmol1, ${script})`
    a.click();
}

// Function that sets the data for the 2D categories plot
// Sets PDB positions for later usage and colors for each peptide
function getCategoriesPeptides2D(data) {
    let traces = [];
    let categories = [];
    let colors = createColors(resultJS.aminoAcidCodes); // Generates colors for each amino acid
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].peptide) === -1) {
            traces.push({
                type: "scatter",
                x: [],
                y: [],
                freetext: [],
                mode: "markers",
                marker: {size: 10, color: []},
                text: `Peptide: ${data[i].peptide}`,
                name: data[i].peptide,
            });
            categories.push(data[i].peptide);
        } 
        traces[categories.indexOf(data[i].peptide)].x.push(data[i].x);
        traces[categories.indexOf(data[i].peptide)].y.push(data[i].y);
        traces[categories.indexOf(data[i].peptide)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
        let rgb = colors.filter(c => c.x === data[i].peptide)[0].y;
        traces[categories.indexOf(data[i].peptide)].marker.color.push(`rgb(${Math.floor(rgb[0])}, ${Math.floor(rgb[1])}, ${Math.floor(rgb[2])})`);
    }
    return traces
};

// Function that sets the data for the 3D categories plot
// Sets PDB positions for later usage and colors for each peptide
function getCategoriesPeptides3D(data) {
    let traces = [];
    let categories = [];
    let colors = createColors(resultJS.aminoAcidCodes);
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].peptide) === -1) {
            traces.push({
                type: "scatter3d",
                x: [],
                y: [],
                z: [],
                freetext: [],
                mode: "markers",
                marker: {size: 2, color: []},
                text: `Peptide: ${data[i].peptide}`,
                name: data[i].peptide,
            });
            categories.push(data[i].peptide);

        }
        traces[categories.indexOf(data[i].peptide)].x.push(data[i].x);
        traces[categories.indexOf(data[i].peptide)].y.push(data[i].y);
        traces[categories.indexOf(data[i].peptide)].z.push(data[i].z);
        traces[categories.indexOf(data[i].peptide)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
        let rgb = colors.filter(c => c.x === data[i].peptide)[0].y;
        traces[categories.indexOf(data[i].peptide)].marker.color.push(`rgb(${Math.floor(rgb[0])}, ${Math.floor(rgb[1])}, ${Math.floor(rgb[2])})`);
    }
    return traces
};

// Function that sets the data for the 2D categories plot
// Sets PDB positions for later usage and colors for each chain
function getCategoriesChains2D(data, colors) {
    let traces = [];
    let categories = [];
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].chain) === -1) {
            traces.push({
                type: "scatter",
                x: [],
                y: [],
                freetext: [],
                mode: "markers",
                marker: {size: 10,
                color: [],
                opacity: 1},
                text: `Chain: ${data[i].chain}`,
                name: data[i].chain,
                visible: false,
            });
            categories.push(data[i].chain);

        } 
        traces[categories.indexOf(data[i].chain)].x.push(data[i].x);
        traces[categories.indexOf(data[i].chain)].y.push(data[i].y);
        traces[categories.indexOf(data[i].chain)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
        let rgb = colors.filter(c => c.x === data[i].chain)[0].y;
        traces[categories.indexOf(data[i].chain)].marker.color.push(`rgb(${Math.floor(rgb[0])}, ${Math.floor(rgb[1])}, ${Math.floor(rgb[2])})`);
    }
    return traces;
};

// Function that sets the data for the 3D categories plot
// Sets PDB positions for later usage and colors for each chain
function getCategoriesChains3D(data, colors) {
    let traces = [];
    let categories = [];
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].chain) === -1) {
            traces.push({
                type: "scatter3d",
                x: [],
                y: [],
                z: [],
                freetext: [],
                mode: "markers",
                marker: {size: 2,
                    opacity: 1,
                color: []
                },
                text: `Chain: ${data[i].chain}`,
                name: data[i].chain,
                visible: false,
            });
            categories.push(data[i].chain);
        } 
        traces[categories.indexOf(data[i].chain)].x.push(data[i].x);
        traces[categories.indexOf(data[i].chain)].y.push(data[i].y);
        traces[categories.indexOf(data[i].chain)].z.push(data[i].z);
        traces[categories.indexOf(data[i].chain)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
        let rgb = colors.filter(c => c.x === data[i].chain)[0].y;
        traces[categories.indexOf(data[i].chain)].marker.color.push(`rgb(${Math.floor(rgb[0])}, ${Math.floor(rgb[1])}, ${Math.floor(rgb[2])})`);
    }
    return traces;
};

// Function that sets the data for the 2D categories plot
// Sets PDB positions for later usage and colors for each secundary structure
function getCategoriesStructure2D(data) {
    let traces = [];
    let categories = [];
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].structure) === -1) {
            traces.push({
                type: "scatter",
                x: [],
                y: [],
                freetext: [],
                mode: "markers",
                marker: {size: 10},
                text: `Structure: ${data[i].structure}`,
                name: data[i].structure,
                visible: false,
            });
            categories.push(data[i].structure);
        }
        traces[categories.indexOf(data[i].structure)].x.push(data[i].x);
        traces[categories.indexOf(data[i].structure)].y.push(data[i].y);
        traces[categories.indexOf(data[i].structure)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
    }
    return traces;
}

// Function that sets the data for the 3D categories plot
// Sets PDB positions for later usage and colors for each secundary structure
function getCategoriesStructure3D(data) {
    let traces = [];
    let categories = [];
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].structure) === -1) {
            traces.push({
                type: "scatter3d",
                x: [],
                y: [],
                z: [],
                freetext: [],
                mode: "markers",
                marker: {size: 2},
                text: `Structure: ${data[i].structure}`,
                name: data[i].structure,
                visible: false,
            });
            categories.push(data[i].structure);
        }
        traces[categories.indexOf(data[i].structure)].x.push(data[i].x);
        traces[categories.indexOf(data[i].structure)].y.push(data[i].y);
        traces[categories.indexOf(data[i].structure)].z.push(data[i].z);
        traces[categories.indexOf(data[i].structure)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
    }
    return traces;
}

function getCategoriesStandard2D(data) {
    return [{
        type: "scatter",
        x: data.x,
        y: data.y,
        mode: "markers",
        marker: {size: 10, color: 'rgb(110, 110, 110)', opacity: 0.5},
        text: `Standard`,
        name: `Standard`
    }];
}

function getCategoriesStandard3D(data) {
    return [{
        type: "scatter3d",
        x: data.x,
        y: data.y,
        z: data.z,
        mode: "markers",
        marker: {size: 2, color: 'rgb(110, 110, 110)', opacity: 0.4},
        text: `Standard`,
        name: `Standard`
    }];
}

// Gets all data for the 2D plot
function getDataPCA2D(json, colors, scores) {
    let tracesPeptides = getCategoriesPeptides2D(json, colors);
    let tracesChains = getCategoriesChains2D(json, colors);
    let tracesStructure = getCategoriesStructure2D(json);
    let standard = getCategoriesStandard2D(scores);


    let buttonVisiblePeptide2D = Array(tracesPeptides.length).fill(true).concat(Array(tracesChains.length).fill(false)).concat(Array(tracesStructure.length).fill(false), true);
    let buttonVisibleChain2D = Array(tracesPeptides.length).fill(false).concat(Array(tracesChains.length).fill(true)).concat(Array(tracesStructure.length).fill(false), true);
    let buttonVisibleStrucutre2D = Array(tracesPeptides.length).fill(false).concat(Array(tracesChains.length).fill(false)).concat(Array(tracesStructure.length).fill(true), true);

    return {
        "data": tracesPeptides.concat(tracesChains, tracesStructure, standard).reverse(),
        "initialView": buttonVisiblePeptide2D.reverse(),
        "secondaryView": buttonVisibleChain2D.reverse(),
        "thirdView": buttonVisibleStrucutre2D.reverse()
    }
};

// Gets all data for the 3D plot
function getDataPCA3D(json, colors, scores) {
    let tracesPeptides = getCategoriesPeptides3D(json, colors);
    let tracesChains = getCategoriesChains3D(json, colors);
    let tracesStructure = getCategoriesStructure3D(json, colors);
    let standard = getCategoriesStandard3D(scores);

    let buttonVisiblePeptide3D = Array(tracesPeptides.length).fill(true).concat(Array(tracesChains.length).fill(false)).concat(Array(tracesStructure.length).fill(false), true);
    let buttonVisibleChain3D = Array(tracesPeptides.length).fill(false).concat(Array(tracesChains.length).fill(true)).concat(Array(tracesStructure.length).fill(false), true);
    let buttonVisibleStrucutre3D = Array(tracesPeptides.length).fill(false).concat(Array(tracesChains.length).fill(false)).concat(Array(tracesStructure.length).fill(true), true);

    return {
        "data": tracesPeptides.concat(tracesChains, tracesStructure, standard),
        "initialView": buttonVisiblePeptide3D,
        "secondaryView": buttonVisibleChain3D,
        "thirdView": buttonVisibleStrucutre3D
    }
};

// Initialized the plotly buttons
function initializePlotlyButtons(initialView, secondaryView, thirdView) {
    return [{
        buttons: [
            {
                args: [
                    {"visible": initialView},
                    {"title": "Peptides", "showlegend": initialView}
                ],
                label: "Peptides",
                method: "update"
            },
            {
                args: [
                    {"visible": secondaryView},
                    {"title": "Chains", "showlegend": secondaryView}
                ],
                label: "Chains",
                method: "update"
            },
            {
                args: [
                    {"visible": thirdView},
                    {"title": "Structures", "showlegend": thirdView}
                ],
                label: "Structures",
                method: "update"
            }
        ],
        direction: "left",
        pad: {"r": 10, "t": 10},
        showactive: true,
        type: "buttons",
        x: 0.1,
        xanchor: "left",
        y: 1.1,
        yanchor: "top"
    }]
};

// Default layout for the 2D plot
function setLayoutPlotly2D(updatemenus) {
    return  {
        autosize: true,
        title: "Peptides",
        margin: {
            l: 0,
            r: 0,
            b: 0,
            t: 0,
            pad: 4
        },
        hovermode: "closest",
        updatemenus: updatemenus,
        showlegend: true,
        paper_bgcolor:"white",
        plot_bgcolor:"#FFFFFF",
        legend: {
            y: 0.5,
            yref: 'paper',
            font: {
                size: 10,
            },
        }
    }
};

// Default layout for the 3D plot
function setLayoutPlotly3D(updatemenus) {
    return {
        autosize: true,
        title: "Peptides",
        margin: {
            l: 0,
            r: 0,
            b: 0,
            t: 0,
            pad: 4
        },
        hovermode: "closest",
        updatemenus: updatemenus,
        showlegend: true,
        scene: {
            aspectmode: "data",
            xaxis: {
                showspikes: false,
                backgroundcolor: "#edf3fa",
                showbackground: true
            },
            yaxis: {
                showspikes: false,
                backgroundcolor: "#edf3fa",
                showbackground: true
            },
            zaxis: {
                showspikes: false,
                backgroundcolor: "#edf3fa",
                showbackground: true
            },
        },
        paper_bgcolor:"white",
        plot_bgcolor:"#00FF00",
        legend: {
            y: 0.5,
            yref: 'paper',
            font: {
                size: 10,
            },
        },
    }
};

// Default JSMOL script
function getInfoProtein3d(pdb) {
    return {
        width: 400,
        height: 292,
        debug: false,
        j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
        color: "0xC0C0C0",
        disableJ2SLoadMonitor: true,
        disableInitialConsole: true,
        addSelectionOptions: false,
        serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
        use: "HTML5",
        readyFunction: null,
        script: `load "=${pdb}"; cartoons only; color structure; background white; zoom 100`
    }
};

document.getElementById('zoom-btn').addEventListener('click', () => {
    hideGlobal();
})

document.getElementById('reset').addEventListener('click', () => {
    resetScript();
})

document.getElementById("placeholder-scatter").style.display= 'none';
document.getElementById("placeholder-scatter-3d").style.display= 'none';
document.getElementById("place-text").style.display= 'none';
document.addEventListener('DOMContentLoaded', (event) => {
    // This function will call the function to create a temporary file and handles the response
    (async function getData() {
        const tokenResponse = await fetch("/csrf-token", {method: "GET"});
        const csrfToken = await tokenResponse.text();
        const fetchParameters = { method: 'POST', headers: {'X-CSRF-TOKEN': csrfToken}}
        const response = await fetch("/create_temp_file", fetchParameters);
        if (response.ok) {
            const chainResponse = await fetch("/get_chains", fetchParameters);
            const chainResult = await chainResponse.json();

            let value = document.getElementById("pdb-structure").textContent;
            value = value.slice(value.indexOf(":") + 2, value.length);
            if (chainResult["bytes"].length > 2) {
                let colors = createColors(Object.keys(JSON.parse(chainResult["bytes"])));
                setChain(chainResult, colors);
              
                // // Functionality for the 3D protein plot
                if (value != null) {
                    let Info = getInfoProtein3d(value);
                    $("#protein").html(Jmol.getAppletHtml("jmol1", Info))
                }
                const dataResponse = await fetch("/perform_pca_analysis", fetchParameters);
                let dataResult = await dataResponse.json();
                dataResult = JSON.parse(dataResult["bytes"]);
                let scores = dataResult.scores;
                let {scores: _, ...result} = dataResult;
                if (result["error"] !== undefined) {
                    window.location.href = `/pdb_error?code=${value}&message=${result["error"]}`
                    return;
                }

                create3dPlot(result, colors, scores)
                create2dPlot(result, colors, scores)

                document.getElementById("download-link").onclick = function () {
                    let blob = new Blob([JSON.stringify(result, null, 4)], {type: "text/json"})
                    saveAs(blob, "result.json");
                }

            } else {
                window.location.href = `/pdb_error?code=${value}`
            }
        } else {
            console.log(response)
        }
    })();
    // This function sets the metadata of the pdb on top of the site. Like the paper link and stuff
    (function () {
        let value = document.getElementById("pdb-structure").textContent;
        value = value.slice(value.indexOf(":") + 2, value.length);
        const url = `https://data.rcsb.org/rest/v1/core/entry/${value}`;
        fetch(url, { method: 'GET' })
            .then((response) => {
                if (!response.ok) {
                    toastr.error(`${value} is not a valid PDB code!`);
                    return Promise.reject(response);
                }
                return response.json();
            })
            .then((result) => {
                if (result["citation"][0]["pdbx_database_id_doi"] != null) {
                    document.getElementById("paper-link").classList.add("is-underlined");
                    document.getElementById("paper-link").onclick = function (){window.open(`http://dx.doi.org/${result["citation"][0]["pdbx_database_id_doi"]}`, '_blank').focus();};
                } else {
                    document.getElementById("paper-link").textContent = "Paper to be published";
                    document.getElementById('paper-link').removeAttribute('id');
                }

                document.getElementById("pdb-title").textContent = result["struct"]["title"];
                document.getElementById("pdb-subtitle").textContent = result["struct"]["pdbx_descriptor"];
                document.getElementById("pdb-link").onclick = function (){window.open(`https://www.rcsb.org/structure/${result["entry"]["id"]}`, '_blank').focus();};
                document.getElementById('skeleton-loader').removeAttribute('id');
            })
            .catch((error) => console.log(error));
    })();
});