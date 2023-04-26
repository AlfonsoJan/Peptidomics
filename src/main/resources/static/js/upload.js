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

function getFileName(el) {
    document.getElementById("file-selected").textContent = el.files[0].name;
}

// Function that checks if it is valid pdb code
async function checkPDBCode(value) {
    console.log(value);
    const response = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${value}`, { method: 'GET' })
    return response;
}
// Function that checks if it is valid pdb code because some pdb files are
// big that they don't have a pdb format
async function checkPDBSize(value) {
    const response = fetch(`https://files.rcsb.org/download/${value}.pdb`)
    return response;
}
document.addEventListener('DOMContentLoaded', (event) => {
    // Function that checks every thing before submitting the pdb code form
    document.getElementById('pdb-input-form').addEventListener('submit', async function (event) {
        event.preventDefault();
        document.getElementById("pdb-input-form").lastElementChild.classList.add("is-loading");

        // Oligo length check > 1
        if (document.getElementById("oligo-param").value < 2) {
            toastr.error(`Please type in a length higher then 1!`);
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
        document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
        window.location.replace(document.referrer);
        document.getElementById("pdb-input-form").submit();
    });
}, false);