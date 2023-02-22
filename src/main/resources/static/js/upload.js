
document.addEventListener('DOMContentLoaded', (event) => {
    document.querySelectorAll(".switch").forEach((theSwitch) => {
        theSwitch.addEventListener("click", handleClickEvent, false);
    });
    function handleClickEvent(evt) {
        const el = evt.target;

        if (el.getAttribute("aria-checked") === "true") {
            el.setAttribute("aria-checked", "false");
        } else {
            el.setAttribute("aria-checked", "true");
        }
    }
});