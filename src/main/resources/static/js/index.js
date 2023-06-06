// author: Wouter Zeevat

document.getElementById("information-card").addEventListener("click", () => {
    let el = document.getElementById('more_info');
    el.scrollIntoView({behavior: "smooth"});
})

document.getElementById("download-standard").addEventListener("click", (e) => {
    e.preventDefault();
    getData();
})

async function getData() {
    const tokenResponse = await fetch("/csrf-token", {method: "GET"});
    const csrfToken = await tokenResponse.text();
    const fetchParameters = { method: 'GET', headers: {'X-CSRF-TOKEN': csrfToken}}
    const response = await fetch("/api/v1/eigenvectors", fetchParameters)
    if (response.ok) {
        const result = await response.json();
        let text = "length;x;y;z\n";
        result.map(r => {
            let lengthVector = r.length;
            for (let i = 0; i < r.x.length; i++) {
                text += `${lengthVector};${r.x[i]};${r.y[i]};${r.z[i]}\n`;
            }
        })
        let blob = new Blob([text], {type: "text"})
        saveAs(blob, "eigenvectors.txt");
    }
}

