// copy_code.js
document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("pre > code").forEach(codeBlock => {
    /* ---------- build the button ---------- */
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "copy-code-button";
    btn.textContent = "Copy";

    /* ---------- copy logic (+ fallback) ---------- */
    btn.addEventListener("click", async () => {
      const codeText = codeBlock.textContent;
      try {
        await navigator.clipboard.writeText(codeText);
      } catch (err) {
        const ta = document.createElement("textarea");
        ta.value = codeText;
        ta.style.position = "fixed";   // avoid scrolling jump
        document.body.appendChild(ta);
        ta.select();
        document.execCommand("copy");
        document.body.removeChild(ta);
      }
      /* UI feedback */
      btn.textContent = "Copied!";
      setTimeout(() => (btn.textContent = "Copy"), 1500);
    });

    /* ---------- place the button ---------- */
    const pre = codeBlock.parentNode;
    pre.style.position = "relative";
    btn.style.position = "absolute";
    btn.style.top = "4px";
    btn.style.right = "4px";
    pre.appendChild(btn);
  });
});
