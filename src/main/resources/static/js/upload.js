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

function getFileName(el) {
    document.getElementById("file-name").textContent = el.files[0].name;
}
(function() {
    document.getElementById("upload-input").checked = true;
    document.getElementById("pdb-input").checked = false;
    document.getElementById("pdb-input-form").style.display = 'none';
    document.getElementById("file-upload").required = true;
})()
// document.addEventListener('DOMContentLoaded', (event) => {
//     (function resetRadioButtons() {
//         document.getElementById("upload-input").checked = true;
//         document.getElementById("pdb-input").checked = false;
//     })()
// }, false);