// Jan Alfonso
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
function getInputType(el) {
    switch(el.value) {
        case "1":
            document.getElementById("pdb-input-form").style.display = 'none';
            document.getElementById("file-upload-form").style.display = '';
            document.getElementById("input-pdb").required = false;
            document.getElementById("file-upload").required = true;
            break;
        case "2":
            document.getElementById("pdb-input-form").style.display = '';
            document.getElementById("file-upload-form").style.display = 'none';
            document.getElementById("input-pdb").required = true;
            document.getElementById("file-upload").required = false;
            break;
        default:
            break;
    }
}
function getCompare(el) {
    const val = el.value;
    document.getElementById("compareCode").value = val;
    document.getElementById("compareFile").value = val;
}
function getParam(el) {
    const val = el.value;
    document.getElementById("paramCode").value = val;
    document.getElementById("paramFile").value = val;
}
function getFileName(el) {
    document.getElementById("file-name").textContent = el.files[0].name;
}
(function() {
    document.getElementById("upload-input").checked = true;
    document.getElementById("pdb-input").checked = false;
    document.getElementById("pdb-input-form").style.display = 'none';
    document.getElementById("file-upload").required = true;
})()

async function checkPDBCode(value) {
    const response = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${value}`, { method: 'GET' })
    return response;
}
async function checkPDBSize(value) {
    const response = fetch(`https://files.rcsb.org/download/${value}.pdb`)
    return response;
}
document.addEventListener('DOMContentLoaded', (event) => {
    document.getElementById('pdb-input-form').addEventListener('submit', async function (event) {
        event.preventDefault();
        document.getElementById("pdb-input-form").lastElementChild.classList.add("is-loading");
        if (document.getElementById("paramFile").value < 2) {
            toastr.error(`Please type in a length higher then 1!`);
            return;
        }

        let compareCode = document.getElementById("compareCode").value;
        const compareStatus = await checkPDBCode(compareCode)
        if (!compareStatus.ok) {
            toastr.error(`${compareCode} is not a valid PDB code!`);
            document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
            return;
        }
        const compareSize = await checkPDBSize(compareCode)
        if (!compareSize.ok) {
            toastr.warning(`${compareCode} does not exist in a PDB format!`);
            document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
            return;
        }

        let pdbCode = document.getElementById("input-pdb").value;
        const pdbStatus = await checkPDBCode(pdbCode)
        if (!pdbStatus.ok) {
            toastr.error(`${pdbCode} is not a valid PDB code!`);
            document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
            return;
        }
        const pdbSize = await checkPDBSize(pdbCode)
        if (!pdbSize.ok) {
            toastr.warning(`${pdbCode} does not exist in a PDB format!`);
            document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
            return;
        }
        window.location.replace(document.referrer);
        document.getElementById("pdb-input-form").submit();
    });
    document.getElementById('file-upload-form').addEventListener('submit', async function (event) {
        event.preventDefault();
        document.getElementById("file-upload-form").lastElementChild.classList.add("is-loading");
        let length = document.getElementById("paramFile").value;
        if (isNaN(length)) {
            toastr.warning(`${length} is not a number!`);
            document.getElementById("file-upload-form").lastElementChild.classList.remove("is-loading");
            return;
        }
        if (length < 2) {
            toastr.error(`Please type in a length higher then 1!`);
            document.getElementById("file-upload-form").lastElementChild.classList.remove("is-loading");
            return;
        }

        let compareCode = document.getElementById("compareFile").value;
        const compareStatus = await checkPDBCode(compareCode)
        if (!compareStatus.ok) {
            toastr.error(`${compareCode} is not a valid PDB code!`);
            document.getElementById("file-upload-form").lastElementChild.classList.remove("is-loading");
            return;
        }

        let value = document.getElementById("file-upload").value;
        if (value.length > 0) {

            window.location.replace(document.referrer);
            document.getElementById("file-upload-form").submit();
        }
    });
}, false);