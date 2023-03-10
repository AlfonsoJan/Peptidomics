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
document.getElementById("info-pdb").style.display = 'none';
document.addEventListener('DOMContentLoaded', (event) => {
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
                document.getElementById("info-pdb").style.display = '';
                document.getElementById("pdb-title").textContent = result["struct"]["title"];
                document.getElementById("pdb-subtitle").textContent = result["struct"]["pdbx_descriptor"];
                document.getElementById("pdb-link").onclick = function (){window.open(`https://www.rcsb.org/structure/${result["entry"]["id"]}`, '_blank').focus();};
                document.getElementById('skeleton-loader').removeAttribute('id');
            })
            .catch((error) => console.log(error));
    })();
});