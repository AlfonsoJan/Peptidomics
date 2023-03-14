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

                                document.getElementById("placeholder-pca").src = `data:image/png;base64,${pca_result["bytes"]}`;
                                document.getElementById("placeholder-scatter").src = `data:image/png;base64,${scatter_result["bytes"]}`;

                                document.getElementById("placeholder-pca").style.display= '';
                                document.getElementById("placeholder-scatter").style.display= '';
                            })
                    })
            })
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