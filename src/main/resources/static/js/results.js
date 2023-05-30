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
        "UNK", // UNKNOWN
    ],

    points: [],
    hidden: false,
    reset: true,
}

let helperFunctions = {
    // Hides all the things on page load
    setPage() {
        document.getElementById("middle-selector").style.display= 'none';
        document.getElementById("left-selector").style.display= 'none';
        document.getElementById("right-selector").style.display= 'none';
        document.getElementById("placeholder-scatter").style.display= 'none';
        document.getElementById("placeholder-scatter-3d").style.display= 'none';
        document.getElementById("place-text").style.display= 'none';
    },
    // Function that converts hue to RGB
    HSLToRGB (h, s, l) {
        s /= 100;
        l /= 100;
        const k = n => (n + h / 30) % 12;
        const a = s * Math.min(l, 1 - l);
        const f = n =>
            l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
        return [255 * f(0), 255 * f(8), 255 * f(4)];
    },
    // Function that creates perfectly distanced colors that follow the rainbow pattern
    // according to an array of keys, creates {"A": red, "B": green, "C": pink}
    createColors(keys) {
        let colorArray = [];
        for (let i = 0; i < keys.length; i++) {
            let hsl_value = 255 / keys.length * i;
            let rgb = this.HSLToRGB(hsl_value, 100, 70);
            colorArray.push(rgb)
        }

        return keys.map((x, i) => ({x, y: colorArray[i]}));
    },
    // Function that will shuffle an array
    shuffle(array, seed) {                // <-- ADDED ARGUMENT
        let m = array.length, t, i;
        // While there remain elements to shuffle…
        while (m) {
            // Pick a remaining element…
            i = Math.floor(this.random(seed) * m--);        // <-- MODIFIED LINE
            // And swap it with the current element.
            t = array[m];
            array[m] = array[i];
            array[i] = t;
            ++seed                                     // <-- ADDED LINE
        }
        return array;
    },
    // Random function for the shuffle
    random(seed) {
        let x = Math.sin(seed++) * 10000;
        return x - Math.floor(x);
    },
    // GET UNIQUE VALUES FROM ARR
    onlyUnique(value, index, array) {
        return array.indexOf(value) === index;
    },
    // Function that set the chains on the site
    setChain(chains, colors) {
        // Function for the click action
        let clickChain = function (item) {
            let ele = item.target;
            while (!ele.parentElement.classList.contains("card")) {
                ele = ele.parentElement;
            }
            let chain = ele.children[0].textContent.charAt(ele.children[0].textContent.length - 1);
            let color = ele.children[0].style.backgroundColor.replace("rgb(", "").replace(")", "").replace(" ", "");
            jMOLHelpers.resetScript();

            let a = document.createElement("a");
            let script = `"select all; color [90,90,90]; select chain=${chain}; cartoons only; color [${color}]; background white; zoom 0"`;
            a.href = `javascript:Jmol.script(jmol1, ${script})`
            a.click();
            resultJS.reset = true;
        }
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
}

let jMOLHelpers = {
    // Selects a certain part of the JSMOL
    selectView(datapoints) {

        // Unhides the JMOL to rehide it after! BUG FIX
        let was_hidden = false;
        if (resultJS.hidden) {
            this.hideGlobal();
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
            this.hideGlobal();
            script = `"zoom 0"`;
        }
        resultJS.reset = false;
    },
    // Default JSMOL script
    getInfoProtein3d(pdb) {
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
    },
    // Zooms in on the selected part and hides the other parts in JSMOL
    hideGlobal() {
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
            let scriptParts = []
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
    },
    // Resets the JSMOL view
    resetScript() {
        // Unhides if it was hidden to prevent bugs
        if (resultJS.hidden) {
            this.hideGlobal();
        }
        document.getElementById("zoom-btn").classList.add("is-hidden");
        let a = document.createElement("a");
        let script = `"select all; cartoons only; color structure; background white; zoom 0"`;
        a.href = `javascript:Jmol.script(jmol1, ${script})`
        a.click();
        document.getElementById("3d-text").innerText = "0 selected";
        resultJS.reset = true;
    }
}


document.getElementById('zoom-btn').addEventListener('click', () => {
    jMOLHelpers.hideGlobal();
})

document.getElementById('reset').addEventListener('click', () => {
    jMOLHelpers.resetScript();
})



document.getElementById('right').addEventListener('change', function() {
    let right = this.value;
    let middle = document.getElementById('middle').value;
    let left = document.getElementById('left').value;
    PlotContainer.updatePlots(left, middle, right);
})

document.getElementById('middle').addEventListener('change', function() {
    let right = document.getElementById('right').value;
    let middle = this.value;
    let left = document.getElementById('left').value;
    PlotContainer.updatePlots(left, middle, right);
})

document.getElementById('left').addEventListener('change', function() {
    let right = document.getElementById('right').value;
    let middle = document.getElementById('middle').value;
    let left = this.value;
    PlotContainer.updatePlots(left, middle, right);
})

let PlotContainer = {
    seed: 3,
    colorsPeptides: helperFunctions.createColors(resultJS.aminoAcidCodes),
    div2D: "placeholder-scatter",
    div3D: "placeholder-scatter-3d",
    config: {responsive: true},
    set standardData(data) {
        this.standard = data
    },
    set colorArr(colArr) {
        this.colors = colArr;
    },
    initializePlotlyButtons(initialView, secondaryView, thirdView) {
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
    },
    getStandardTraces2D() {
        return [{
            type: "scatter",
            x: this.standard.x,
            y: this.standard.y,
            z: this.standard.z,
            mode: "markers",
            marker: {size: 10, color: 'rgb(110, 110, 110)', opacity: 0.4},
            text: `Standard`,
            name: `Standard`
        }];
    },
    getStandardTraces3D() {
        return [{
            type: "scatter3d",
            x: this.standard.x,
            y: this.standard.y,
            z: this.standard.z,
            mode: "markers",
            marker: {size: 2, color: 'rgb(110, 110, 110)', opacity: 0.4},
            text: `Standard`,
            name: `Standard`
        }];
    },
    set dataPlot(data) {
        this.data = data
    },
    getMetadata() {
        return Object.keys(this.data).map(key => {
            let data = this.data[key]
            return {
                "atomnos": [data.atomnos.min, data.atomnos.max],
                "chains": data.chain,
                "peptide": data.peptide,
                "structure": data.structure
            }
        });
    },
    setUniqueCat() {
        let selectRight = document.getElementById('right');
        let selectMiddle = document.getElementById('middle');
        let selectLeft = document.getElementById('left');
        let uniqueChains = [].concat.apply([], Object.keys(this.data).map(key => {
            return this.data[key].chain
        })).filter(helperFunctions.onlyUnique);
        let uniquePeptides = [].concat.apply([], Object.keys(this.data).map(key => {
            return this.data[key].peptide
        })).filter(helperFunctions.onlyUnique);
        let uniqueStructure= [].concat.apply([], Object.keys(this.data).map(key => {
            return this.data[key].structure
        })).filter(helperFunctions.onlyUnique);
        let uniqueDataList = uniqueChains.concat(uniquePeptides, uniqueStructure);
        let amount = this.data[0].chain.length;
        if (amount === 1) {
            document.getElementById("middle-selector").style.display= '';
        } else if (amount === 2) {
            document.getElementById("left-selector").style.display= '';
            document.getElementById("right-selector").style.display= '';
        } else {
            document.getElementById("middle-selector").style.display= '';
            document.getElementById("left-selector").style.display= '';
            document.getElementById("right-selector").style.display= '';
        }
        uniqueDataList.forEach(element => {
            let optLeft = document.createElement('option');
            optLeft.value = element;
            optLeft.innerHTML = element;
            selectLeft.appendChild(optLeft);

            let optMiddle = document.createElement('option');
            optMiddle.value = element;
            optMiddle.innerHTML = element;
            selectMiddle.appendChild(optMiddle);

            let optRight = document.createElement('option');
            optRight.value = element;
            optRight.innerHTML = element;
            selectRight.appendChild(optRight);
        });
    },
    getViews(dataPeptides, dataChains, dataStructure) {
        return [
            helperFunctions.shuffle(Array(dataPeptides.length).fill(true).concat(Array(dataChains.length).fill(false)).concat(Array(dataStructure.length).fill(false)), this.seed).concat(true),
            helperFunctions.shuffle(Array(dataPeptides.length).fill(false).concat(Array(dataChains.length).fill(true)).concat(Array(dataStructure.length).fill(false)), this.seed).concat(true),
            helperFunctions.shuffle(Array(dataPeptides.length).fill(false).concat(Array(dataChains.length).fill(false)).concat(Array(dataStructure.length).fill(true)), this.seed).concat(true)
        ]

    },
    removeSpinners() {
        let elem3D = document.getElementById("spinner-scatter-3d");
        elem3D.parentNode.removeChild(elem3D);
        document.getElementById(this.div3D).style.display= '';

        let elem2D = document.getElementById("spinner-scatter");
        elem2D.parentNode.removeChild(elem2D);
        document.getElementById(this.div2D).style.display= '';
    },
    setInitialPlots() {
        let [tracesPeptides2D, tracesChains2D, tracesStructure2D, tracesPeptides3D, tracesChains3D, tracesStructure3D, standard2D, standard3D] = this.getAllData(this.data);
        // 2D and 3D shuffled data
        let data2D = helperFunctions.shuffle(tracesPeptides2D.concat(tracesChains2D, tracesStructure2D), this.seed).concat(standard2D).reverse();
        let data3D = helperFunctions.shuffle(tracesPeptides3D.concat(tracesChains3D, tracesStructure3D), this.seed).concat(standard3D);
        // 2D VIEWS
        let [initialView2D, secondaryView2D, thirdView2D] = this.getViews(tracesPeptides2D.reverse(), tracesChains2D.reverse(), tracesStructure2D.reverse());
        let [initialView3D, secondaryView3D, thirdView3D] = this.getViews(tracesPeptides3D, tracesChains3D, tracesStructure3D);
        // REMOVE SPINNERS
        this.removeSpinners();
        // BUTTONS
        let updateMenus2D = this.initializePlotlyButtons(initialView2D.reverse(), secondaryView2D.reverse(), thirdView2D.reverse());
        let updateMenus3D = this.initializePlotlyButtons(initialView3D, secondaryView3D, thirdView3D);
        // 2D Plot
        Plotly.newPlot(this.div2D, data2D, this.get2DLayout(updateMenus2D), this.config);
        // 3D PLOT
        Plotly.newPlot(this.div3D, data3D, this.get3DLayout(updateMenus3D), this.config);
        // CLICK FUNCTION
        let myPlot2D = document.getElementById(this.div2D);
        myPlot2D.on("plotly_click", this.clickPoints);
        myPlot2D.on("plotly_selected", this.selectMultiplePoints);
        let myPlot3D = document.getElementById(this.div3D);
        myPlot3D.on("plotly_click", this.clickPoints);
    },
    selectMultiplePoints(datapoints) {
        let points = [];
        datapoints.points.forEach(data => {
            let index = data.pointIndex;
            if (!data.data.hasOwnProperty("freetext")) return; // If it's a gray dot AKA background data
            let point = data.data.freetext[index][0].atomnos
            points.push(point);
        });
        jMOLHelpers.selectView(points);
    },
    clickPoints(datapoints) {
        if (!datapoints.points[0].data.hasOwnProperty("freetext")) return;
        jMOLHelpers.selectView([datapoints.points[0].data.freetext[datapoints.points[0].pointNumber][0].atomnos])
    },
    get3DLayout(updatemenus) {
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
            plot_bgcolor:"#FFFFFF",
            legend: {
                y: 0.5,
                yref: 'paper',
                font: {
                    size: 10,
                },
            }
        }
    },
    get2DLayout(updatemenus) {
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
    },
    updatePlots(left, middle, right) {
        jMOLHelpers.resetScript();
        let metaData = this.getMetadata();
        let newDataIndex = Object.keys(metaData).map(key => {
            return metaData[key]
        }).reduce(function(acc, d, index) {
            if (d.chains.length === 2) {
                let rightIndex = (d.chains[0] === right | d.peptide[0] === right | d.structure[0] === right)
                let leftIndex = (d.chains[1] === left | d.peptide[1] === left | d.structure[1] === left)
                if (right === "_") rightIndex = true;
                if (left === "_") leftIndex = true;
                if (leftIndex && rightIndex) acc.push(index);
            } else if (d.chains.length === 3) {
                let rightIndex = (d.chains[0] === right | d.peptide[0] === right | d.structure[0] === right)
                let middleIndex = (d.chains[1] === middle | d.peptide[1] === middle | d.structure[1] === middle)
                let leftIndex = (d.chains[2] === left | d.peptide[2] === left | d.structure[2] === left)
                if (right === "_") rightIndex = true;
                if (middle === "_") middleIndex = true;
                if (left === "_") leftIndex = true;
                if (leftIndex && rightIndex && middleIndex) acc.push(index);
            } else {
                let middleIndex = (d.chains[0] === middle | d.peptide[0] === middle | d.structure[0] === middle);
                if (middle === "_") middleIndex = true;
                if (middleIndex) acc.push(index);
            }
            return acc;
        }, []);
        let newJson = newDataIndex.map(key => this.data[key])
        if (newJson.length < 1) {
            toastr.info("This combination does not exist in this plot!");
            return;
        }
        let [tracesPeptides2D, tracesChains2D, tracesStructure2D, tracesPeptides3D, tracesChains3D, tracesStructure3D, standard2D, standard3D] = this.getAllData(newJson);
        // 2D and 3D shuffled data
        let data2D = helperFunctions.shuffle(tracesPeptides2D.concat(tracesChains2D, tracesStructure2D), this.seed).concat(standard2D).reverse();
        let data3D = helperFunctions.shuffle(tracesPeptides3D.concat(tracesChains3D, tracesStructure3D), this.seed).concat(standard3D);
        // 2D and 3D VIEWS
        let [initialView2D, secondaryView2D, thirdView2D] = this.getViews(tracesPeptides2D, tracesChains2D, tracesStructure2D);
        let [initialView3D, secondaryView3D, thirdView3D] = this.getViews(tracesPeptides3D, tracesChains3D, tracesStructure3D);
        // BUTTONS
        let updateMenus2D = this.initializePlotlyButtons(initialView2D.reverse(), secondaryView2D.reverse(), thirdView2D.reverse());
        let updateMenus3D = this.initializePlotlyButtons(initialView3D, secondaryView3D, thirdView3D);
        // 2D Plot
        Plotly.react(this.div2D, data2D, this.get2DLayout(updateMenus2D), this.config);
        // 3D PLOT
        Plotly.react(this.div3D, data3D, this.get3DLayout(updateMenus3D), this.config);
        // CLICK FUNCTION
        let myPlot2D = document.getElementById(this.div2D);
        myPlot2D.on("plotly_click", this.clickPoints)
        myPlot2D.on("plotly_selected", this.selectMultiplePoints);
        let myPlot3D = document.getElementById(this.div3D);
        myPlot3D.on("plotly_click", this.clickPoints)
    },
    getAllData(data) {
        // 2D DATA
        let tracesPeptides2D = this.getCategoriesDATA(data, "peptide", true, "2D");
        let tracesChains2D = this.getCategoriesDATA(data, "chain", false, "2D");
        let tracesStructure2D = this.getCategoriesDATA(data, "structure", false, "2D");
        // 3D DATA
        let tracesPeptides3D = this.getCategoriesDATA(data, "peptide", true, "3D");
        let tracesChains3D = this.getCategoriesDATA(data, "chain", false, "3D")
        let tracesStructure3D = this.getCategoriesDATA(data, "structure", false, "3D");
        // STANDARD
        let standard2D = this.getStandardTraces2D();
        let standard3D = this.getStandardTraces3D();
        return [tracesPeptides2D, tracesChains2D, tracesStructure2D, tracesPeptides3D, tracesChains3D, tracesStructure3D, standard2D, standard3D];
    },
    getCategoriesDATA(data, select, visible, view) {
        let metaData = Object.keys(data).map(key => {
            return {
                "atomnos": [data[key].atomnos.min, data[key].atomnos.max],
                "chains": data[key].chain,
                "peptide": data[key].peptide,
                "structure": data[key].structure
            }
        })
        let traces = [];
        let categories = [];
        for (let i = 0; i < Object.keys(data).length; i ++) {
            let cat = data[i][select];

            if (cat.length === 2) {
                cat = cat[0]
            } else if (cat.length === 3) {
                cat = cat[1]
            }
            if (categories.indexOf(cat) === -1) {
                traces.push({
                    x: [],
                    y: [],
                    z: [],
                    freetext: [],
                    mode: "markers",
                    marker: {size: 2},
                    text: `${select}: ${cat}`,
                    name: cat,
                })
                categories.push(cat);
                if (!(visible)) {
                    traces[traces.length - 1].visible = false
                }
                if (view === "3D") {
                    traces[traces.length - 1].type = "scatter3d"
                    traces[traces.length - 1].marker.size = 2
                } else {
                    traces[traces.length - 1].type = "scatter"
                    traces[traces.length - 1].marker.size = 10
                }
                if (select === "chain" || select === "peptide") traces[traces.length - 1].marker.color = []
            }

            traces[categories.indexOf(cat)].x.push(data[i].x);
            traces[categories.indexOf(cat)].y.push(data[i].y);
            traces[categories.indexOf(cat)].z.push(data[i].z);
            if (select === "chain") {
                if (cat !== "R") {
                    let rgb = this.colors.filter(c => c.x === cat)[0].y;
                    traces[categories.indexOf(cat)].marker.color.push(`rgb(${Math.floor(rgb[0])}, ${Math.floor(rgb[1])}, ${Math.floor(rgb[2])})`);
                } else {
                    traces[categories.indexOf(cat)].marker.color.push(`rgb(70, 70, 70)`);
                }
            } else if (select === "peptide"){
                let rgb = this.colorsPeptides.filter(c => c.x === cat)[0].y
                traces[categories.indexOf(cat)].marker.color.push(`rgb(${Math.floor(rgb[0])}, ${Math.floor(rgb[1])}, ${Math.floor(rgb[2])})`);
            }

            traces[categories.indexOf(cat)].freetext.push([metaData[i]]);
        }
        return traces;
    },
}

helperFunctions.setPage();

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
                let colors = helperFunctions.createColors(Object.keys(JSON.parse(chainResult["bytes"])));
                helperFunctions.setChain(chainResult, colors);
              
                // Functionality for the 3D protein plot
                if (value != null) {
                    let info = jMOLHelpers.getInfoProtein3d(value);
                    $("#protein").html(Jmol.getAppletHtml("jmol1", info))
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
                PlotContainer.dataPlot = result;
                PlotContainer.standardData = scores;
                PlotContainer.colorArr = colors;
                PlotContainer.setUniqueCat();
                PlotContainer.setInitialPlots();


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