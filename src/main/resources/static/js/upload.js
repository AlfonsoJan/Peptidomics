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
function getParam(el) {
    const val = el.value;
    document.getElementById("param_code").value = val;
    document.getElementById("param_file").value = val;
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
document.addEventListener('DOMContentLoaded', (event) => {
    document.getElementById('pdb-input-form').addEventListener('submit', function(event){
        event.preventDefault();
        document.getElementById("pdb-input-form").lastElementChild.classList.add("is-loading");
        let value = document.getElementById("input-pdb").value;
        fetch(`https://data.rcsb.org/rest/v1/core/entry/${value}`, { method: 'GET' })
            .then((response) => {
                if (!response.ok) {
                    document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
                    toastr.error(`${value} is not a valid PDB code!`);
                    return Promise.reject(response);
                }
                return response.json();
            })
            .then((result) => {
                fetch(`https://files.rcsb.org/download/${value}.pdb`).then(res => {
                    if (!res.ok) {
                        document.getElementById("pdb-input-form").lastElementChild.classList.remove("is-loading");
                        toastr.warning(`${value} does not exist in a PDB format!`);
                        return Promise.reject(res);
                    }
                    document.getElementById("input-pdb").value = result["entry"].id;
                    window.location.replace(document.referrer);
                    document.getElementById("pdb-input-form").submit();
                })
            })
            .catch((error) => console.log(error));
    });
    document.getElementById('file-upload-form').addEventListener('submit', function(event){
        event.preventDefault();
        let value = document.getElementById("file-upload").value;
        if (value.length > 0) {
            document.getElementById("file-upload-form").lastElementChild.classList.add("is-loading");
            window.location.replace(document.referrer);
            document.getElementById("file-upload-form").submit();
        }
    });
}, false);