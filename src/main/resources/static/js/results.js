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

// Function that converts hue to RGB
const HSLToRGB = (h, s, l) => {
    s /= 100;
    l /= 100;
    const k = n => (n + h / 30) % 12;
    const a = s * Math.min(l, 1 - l);
    const f = n =>
        l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
    return [255 * f(0), 255 * f(8), 255 * f(4)];
};

function rgbToHex(r, g, b) {
    return '#' + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
}

let colorArray = [];
for (let i = 0; i < 26; i++) {
    let hsl_value = 255 / 26 * i;
    let rgb = HSLToRGB(hsl_value, 100, 80);
    let hex = rgbToHex(rgb[0], rgb[1], rgb[2])
    colorArray.push(hex)
}
for (let i = colorArray.length - 1; i > 0; i--) {
    let j = Math.floor(Math.random() * (i + 1));
    [colorArray[i], colorArray[j]] = [colorArray[j], colorArray[i]];
}

Plotly.setPlotConfig({
    colors: colorArray
});

let letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split('')
chain_dict = letters.map((x, i) => ({ x, y: colorArray[i] }));

// Function that create the 2d scatter plotly plot for the dimension
function createDimPlot(result) {
    let elem = document.getElementById("spinner-pca");
    elem.parentNode.removeChild(elem);
    document.getElementById("placeholder-pca").style.display= '';
    let trace1 = {
        type: 'scatter',
        x: JSON.parse(result["bytes"])["dim"].x,
        y: JSON.parse(result["bytes"])["dim"].y,
        mode: 'markers',
        marker: {
            color: 'rgb(17, 157, 255)',
            size: 10
        }
    };
    let data = [ trace1 ];
    let layout = {
        autosize: true,
        margin: {
            l: 0,
            r: 0,
            b: 0,
            t: 0,
            pad: 4

        },
        xaxis: {
            range: [0, 10],
        },
        paper_bgcolor:"white",
        plot_bgcolor:"#FFFFFF",
    };
    let config = {responsive: true}
    Plotly.newPlot('placeholder-pca',data,layout,config);
}
// Function that create the 3d scatter plotly plot
function create3dPlot(result) {
    let elem = document.getElementById("spinner-scatter-3d");
    elem.parentNode.removeChild(elem);
    document.getElementById("placeholder-scatter-3d").style.display= '';
    let { data, initialView, secondaryView } = getDataPCA3D(result);
    let updateMenus = initializePlotlyButtons(initialView, secondaryView);
    let layout = setLayoutPlotly3D(updateMenus);
    let config = {responsive: true};
    Plotly.newPlot("placeholder-scatter-3d", data, layout, config);
    let myPlot = document.getElementById("placeholder-scatter-3d");
    myPlot.on("plotly_click", function(datapoints){
        const data = datapoints.points[0].data.freetext[datapoints.points[0].pointNumber];
        const atomMin = parseInt(data[0]);
        const atomMax = parseInt(data[1]);
        let a = document.createElement("a");
        let script = `"spacefill off; select all; cartoons; color [84,84,84]; select atomno>${atomMin - 1} and atomno<${atomMax + 1}; spacefill; color [10,0,255];"`
        a.href = `javascript:Jmol.script(jmol1, ${script})`
        a.click();
    });
}
// Function that create the 2d scatter plotly plot for the PCA results
function create2dPlot(result) {
    let elem = document.getElementById("spinner-scatter");
    elem.parentNode.removeChild(elem);
    document.getElementById("placeholder-scatter").style.display= '';
    let { data, initialView, secondaryView } = getDataPCA2D(result);
    let updateMenus = initializePlotlyButtons(initialView, secondaryView);
    let layout = setLayoutPlotly2D(updateMenus);
    let config = {responsive: true};
    Plotly.newPlot('placeholder-scatter', data, layout, config);
    let myPlot = document.getElementById("placeholder-scatter");
    myPlot.on("plotly_click", function(datapoints){
        const data = datapoints.points[0].data.freetext[datapoints.points[0].pointNumber];
        const atomMin = parseInt(data[0]);
        const atomMax = parseInt(data[1]);
        let a = document.createElement("a");
        let script = `"spacefill off; select all; cartoons; color [84,84,84]; select atomno>${atomMin - 1} and atomno<${atomMax + 1}; spacefill; color [10,0,255];"`
        a.href = `javascript:Jmol.script(jmol1, ${script})`
        a.click();
    });
}
// Function that set the chains on the site
function setChain(chains) {
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
            card.className = "card";
            column.appendChild(card);

            const cardContent = document.createElement("div");
            cardContent.className = "chain-content card-content";
            card.appendChild(cardContent);

            // Loops rainbow colors to get best colors that the plot will late be using
            let hex = chain_dict.filter(word => word.x === c[0])[0]["y"];
            const chainId = document.createElement("p");
            chainId.className = "is-size-5 has-text-weight-bold";
            chainId.textContent = `Chain: ${c[0]}`;
            chainId.style.backgroundColor = `${hex}`;
            cardContent.appendChild(chainId);

            const atomLength = document.createElement("p");
            atomLength.className = "is-size-5";
            atomLength.textContent = `ATOM: ${c[1]} residues`;
            cardContent.appendChild(atomLength);
        });
        document.getElementById("stats-pdb").appendChild(columns);
    }
}
const getCategoriesPepties2D = (data) => {
    let traces = [];
    let categories = [];
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].peptide) === -1) {
            traces.push({
                type: "scatter",
                x: [],
                y: [],
                freetext: [],
                mode: "markers",
                marker: {size: 10},
                text: `Peptide: ${data[i].peptide}`,
                name: data[i].peptide,
            });
            categories.push(data[i].peptide);
        } else {
            traces[categories.indexOf(data[i].peptide)].x.push(data[i].x);
            traces[categories.indexOf(data[i].peptide)].y.push(data[i].y);
            traces[categories.indexOf(data[i].peptide)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
        }
    }
    return traces
};
const getCategoriesPepties3D = (data) => {
    let traces = [];
    let categories = [];
    for (let i = 0; i < Object.keys(data).length; i ++) {
        if (categories.indexOf(data[i].peptide) === -1) {
            traces.push({
                type: "scatter3d",
                x: [],
                y: [],
                z: [],
                freetext: [],
                mode: "markers",
                marker: {size: 2},
                text: `Peptide: ${data[i].peptide}`,
                name: data[i].peptide,
            });
            categories.push(data[i].peptide);
        } else {
            traces[categories.indexOf(data[i].peptide)].x.push(data[i].x);
            traces[categories.indexOf(data[i].peptide)].y.push(data[i].y);
            traces[categories.indexOf(data[i].peptide)].z.push(data[i].z);
            traces[categories.indexOf(data[i].peptide)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
        }
    }
    return traces
};

const getCategoriesChains2D = (data) => {
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
                marker: {size: 10},
                text: `Chain: ${data[i].chain}`,
                name: data[i].chain,
                visible: false,
            });
            categories.push(data[i].chain);
        } else {
            traces[categories.indexOf(data[i].chain)].x.push(data[i].x);
            traces[categories.indexOf(data[i].chain)].y.push(data[i].y);
            traces[categories.indexOf(data[i].chain)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
        }
    }
    return traces;
};
const getCategoriesChains3D = (data) => {
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
                marker: {size: 2},
                color: [],
                text: `Chain: ${data[i].chain}`,
                name: data[i].chain,
                visible: false,
            });
            categories.push(data[i].chain);
        } else {
            traces[categories.indexOf(data[i].chain)].x.push(data[i].x);
            traces[categories.indexOf(data[i].chain)].y.push(data[i].y);
            traces[categories.indexOf(data[i].chain)].z.push(data[i].z);
            traces[categories.indexOf(data[i].chain)].freetext.push([data[i].atomnos.min, data[i].atomnos.max]);
            let hex = chain_dict.filter(c => c.x === data[i].chain)[0]['y'];
            traces[categories.indexOf(data[i].chain)].color.push(`${hex}`);
        }
    }
    console.log(traces);
    return traces;
};
const getDataPCA2D = (json) => {
    let tracesPeptides = getCategoriesPepties2D(json);
    let tracesChains = getCategoriesChains2D(json);

    let buttonVisible = Array(tracesPeptides.length).fill(true).concat(Array(tracesChains.length).fill(false));
    let buttonVisibleReverse = Array(tracesPeptides.length).fill(false).concat(Array(tracesChains.length).fill(true));

    return {
        "data": tracesPeptides.concat(tracesChains),
        "initialView": buttonVisible,
        "secondaryView": buttonVisibleReverse,
    }
};
const getDataPCA3D = (json) => {
    let tracesPeptides = getCategoriesPepties3D(json);
    let tracesChains = getCategoriesChains3D(json);

    let buttonVisible = Array(tracesPeptides.length).fill(true).concat(Array(tracesChains.length).fill(false));
    let buttonVisibleReverse = Array(tracesPeptides.length).fill(false).concat(Array(tracesChains.length).fill(true));

    return {
        "data": tracesPeptides.concat(tracesChains),
        "initialView": buttonVisible,
        "secondaryView": buttonVisibleReverse,
    }
};
const initializePlotlyButtons = (initialView, secondaryView) => {
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
const setLayoutPlotly2D = (updatemenus) => {
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
const setLayoutPlotly3D = (updatemenus) => {
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
        colorway: colorArray
    }
};
const getInfoProtein3d = (pdb) => {
    return {
        width: 400,
        height: 400,
        debug: false,
        j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
        color: "0xC0C0C0",
        disableJ2SLoadMonitor: true,
        disableInitialConsole: true,
        addSelectionOptions: false,
        serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
        use: "HTML5",
        readyFunction: null,
        script: `load "=${pdb}"; cartoons only; color structure; zoom 50; wireframe;`
    }
};
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

            setChain(chainResult)

            let value = document.getElementById("pdb-structure").textContent;
            value = value.slice(value.indexOf(":") + 2, value.length);
            // Functionality for the 3D protein plot
            if (value != null) {
                let Info = getInfoProtein3d(value);
                $("#protein").html(Jmol.getAppletHtml("jmol1", Info))
            }
            const dataResponse = await fetch("/perform_pca_analysis", fetchParameters);
            let dataResult = await dataResponse.json();
            dataResult = JSON.parse(dataResult["bytes"]);
            create3dPlot(dataResult)
            create2dPlot(dataResult)
            //createDimPlot(dataResult)
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