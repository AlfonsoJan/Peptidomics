// author: Jan Alfonso Busker
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

// Function that checks if it is valid pdb code
async function checkPDBCode(value) {
    const response = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${value}`, { method: 'GET' })
    return response;
}
// Function that checks if it is valid pdb code because some pdb files are
// big that they don't have a pdb format
async function checkPDBSize(value) {
    const response = fetch(`https://files.rcsb.org/download/${value}.pdb`)
    return response;
}

document.getElementById("selected").addEventListener('change', (element) => {
    document.getElementById("file-selected").textContent = element.target.files[0].name;
});

document.addEventListener('DOMContentLoaded', (event) => {
    // Function that checks every thing before submitting the pdb code form
    document.getElementById('pdb-input-form').addEventListener('submit', async function (event) {
        event.preventDefault();
        document.getElementById("pdb-input-form").lastElementChild.classList.add("is-loading");

        // Oligo length check > 1
        if (document.getElementById("oligoParam").value < 1 || document.getElementById("oligoParam").value > 30) {
            toastr.error(`Please type in a length higher than 0 and lower than 31!`);
            document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
            return;
        }

        // Check if oligo code or file is set <- Code overrides File
        let pdbCode = document.getElementById("input-pdb").value;
        if (pdbCode.length < 1) {

            if (!document.getElementById("file-selected").textContent.endsWith(".pdb")) {
                toastr.warning(`Please fill in either a good code or a valid file!`);
                document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
                return;
            }
            document.getElementById("pdb-input-form").action = "result_from_files";
        } else {
            const pdbStatus = await checkPDBCode(pdbCode)
            const pdbSize = await checkPDBSize(pdbCode)

            if (!pdbStatus.ok) {
                toastr.error(`${pdbCode} is not a valid PDB code!`);
                document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
                return;
            } else if (!pdbSize.ok) {
                toastr.warning(`${pdbCode} does not exist in a PDB format!`);
                document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
                return;
            }
            document.getElementById("pdb-input-form").action = "result_from_code";
        }
        window.location.replace(document.referrer);
        document.getElementById("pdb-input-form").submit();
    });
}, false);