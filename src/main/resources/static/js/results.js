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
function insertAfter(newNode, existingNode) {
    existingNode.parentNode.insertBefore(newNode, existingNode.nextSibling);
}
document.getElementById("placeholder-pca").style.display= 'none';
document.getElementById("placeholder-scatter").style.display= 'none';
document.getElementById("placeholder-scatter-3d").style.display= 'none';
document.addEventListener('DOMContentLoaded', (event) => {
    (function () {
        fetch("/create_temp_file", { method: 'POST' })
            .then()
            .then(() => {
                fetch("/create_pca_plot", { method: 'POST' })
                    .then(res => res.json())
                    .then(pca_result => {
                        fetch("/create_scatter_plot", { method: 'POST' })
                            .then(res => res.json())
                            .then(scatter_result => {
                                let elem = document.getElementById("spinner-pca");
                                elem.parentNode.removeChild(elem);
                                elem = document.getElementById("spinner-scatter");
                                elem.parentNode.removeChild(elem);

                                document.getElementById("placeholder-pca").style.display= '';
                                document.getElementById("placeholder-scatter").style.display= '';
                                console.log(scatter_result)
                                let trace1 = {
                                    type: 'scatter',
                                    x: JSON.parse(pca_result["bytes"]).x,
                                    y: JSON.parse(pca_result["bytes"]).y,
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
                                        range: [0, 10]  // to set the xaxis range to 0 to 1
                                    },
                                    paper_bgcolor:"white",
                                    plot_bgcolor:"#FFFFFF"
                                };
                                let config = {responsive: true}
                                Plotly.newPlot('placeholder-pca',data,layout,config);


                                trace1 = {
                                    type: 'scatter',
                                    x: JSON.parse(scatter_result["bytes"]).x,
                                    y: JSON.parse(scatter_result["bytes"]).y,
                                    mode: 'markers',
                                    marker: {
                                        color: 'rgba(17, 157, 255, 0.35)',
                                        size: 10
                                    }
                                };
                                data = [ trace1 ];
                                layout = {
                                    autosize: true,
                                    margin: {
                                        l: 0,
                                        r: 0,
                                        b: 0,
                                        t: 0,
                                        pad: 4
                                    },
                                    paper_bgcolor:"white",
                                    plot_bgcolor:"#FFFFFF"
                                };
                                config = {responsive: true}
                                Plotly.newPlot('placeholder-scatter',data,layout,config);
                            })
                    })
                fetch("/pca_plotly_plot", { method: 'POST' })
                    .then(res => res.json())
                    .then(result => {
                        let elem = document.getElementById("spinner-scatter-3d");
                        elem.parentNode.removeChild(elem);
                        document.getElementById("placeholder-scatter-3d").style.display= '';
                        let trace1 = {
                            type: 'scatter3d',
                            x: JSON.parse(result["bytes"]).x,
                            y: JSON.parse(result["bytes"]).y,
                            z: JSON.parse(result["bytes"]).z,
                            mode: 'markers',
                            marker: {
                                size: 2
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
                                }
                            },
                            paper_bgcolor:"white",
                            plot_bgcolor:"#00FF00"
                        };
                        let config = {responsive: true}
                        Plotly.newPlot('placeholder-scatter-3d',data,layout,config);
                    })
            })
    })();
    (function () {
        fetch("/get_stats_pdb", { method: 'POST' })
            .then(response => response.json())
            .then((res) => {
                let loader = document.getElementById("loader-stats");
                loader.parentNode.removeChild(loader);

                let elem = document.getElementById("pdb-stats");
                res["chainList"].forEach(function (i) {
                    const column = document.createElement("div");
                    column.className = "column";

                    const card = document.createElement("div");
                    card.className = "card";
                    column.appendChild(card);

                    const cardContent = document.createElement("div");
                    cardContent.className = "card-content";
                    card.appendChild(cardContent);

                    const chainId = document.createElement("p");
                    chainId.className = "is-size-5 has-text-weight-bold";
                    chainId.textContent = `Chain: ${i["chainId"]}`;
                    cardContent.appendChild(chainId);

                    const seqresLength = document.createElement("p");
                    seqresLength.className = "is-size-5";
                    seqresLength.textContent = `SEQRES: ${i["seqres"].length}`;
                    cardContent.appendChild(seqresLength);

                    const atomLength = document.createElement("p");
                    atomLength.className = "is-size-5";
                    atomLength.textContent = `ATOM: ${i["count"]} residues`;
                    cardContent.appendChild(atomLength);

                    elem.appendChild(column);
                });
            })
            .catch((error) => console.log(error));
    })();
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