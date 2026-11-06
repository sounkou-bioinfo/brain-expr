"""Lance (via le protocole CDP) l'app Electron construite par shinyelectron, attend que webR démarre
l'app Shiny, parcourt les onglets et prend des captures.
Usage : xvfb-run -a python3 tests/test_electron.py <linux-unpacked/exe> <dossier_captures>"""
import sys, time, os, subprocess
from playwright.sync_api import sync_playwright

exe, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
logs = []
with sync_playwright() as p:
    proc = subprocess.Popen([exe, "--no-sandbox", "--remote-debugging-port=9333"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    browser = None
    for _ in range(60):
        try:
            browser = p.chromium.connect_over_cdp("http://127.0.0.1:9333"); break
        except Exception:
            time.sleep(1)
    while not browser.contexts or not browser.contexts[0].pages:
        time.sleep(1)
    win = browser.contexts[0].pages[0]
    win.on("console", lambda m: logs.append(f"[{m.type}] {m.text}"))
    win.set_viewport_size({"width": 1280, "height": 860})
    t0 = time.time()
    frame = None
    # shinylive rend l'app dans une iframe ; on attend le titre de l'app
    while time.time() - t0 < 300:
        for f in win.frames:
            try:
                if f.locator("text=Expression de gènes").count() and f.locator("#dot img, #dot canvas").count():
                    frame = f
                    break
            except Exception:
                pass
        if frame:
            break
        time.sleep(3)
    print(f"app prête : {frame is not None} en {time.time()-t0:.0f}s")
    win.screenshot(path=f"{out}/1_dotplot.png")
    if frame:
        for tab, name in [("Coexpression", "2_coex"), ("Dans le temps", "3_temps"), ("Table", "4_table")]:
            frame.get_by_role("tab", name=tab).click()
            time.sleep(6)
            win.screenshot(path=f"{out}/{name}.png")
    browser.close(); proc.terminate()
print("\n".join([l for l in logs if "error" in l.lower()][:20]))
