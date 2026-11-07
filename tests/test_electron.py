"""Launches the built desktop app, waits for webR to start the Shiny app, opens
every tab and fails if a plot does not render or an output shows an error.

Usage: python tests/test_electron.py <unpacked_dir> <shots_dir>
  <unpacked_dir>: build/electron-app/dist/{linux-unpacked,win-unpacked,mac*}
On Linux without a display, run it under xvfb-run.
"""
import glob, os, subprocess, sys, time
from playwright.sync_api import sync_playwright

unpacked, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)


def find_exe(d):
    if sys.platform == "darwin":
        hits = glob.glob(os.path.join(d, "**", "*.app", "Contents", "MacOS", "*"), recursive=True)
    elif sys.platform == "win32":
        hits = [p for p in glob.glob(os.path.join(d, "*.exe")) if "uninstall" not in p.lower()]
    else:
        hits = [os.path.join(d, "brain-gene-expression")]
    hits = [h for h in hits if os.path.isfile(h)]
    if not hits:
        raise SystemExit(f"no executable found under {d}")
    return hits[0]


def app_frame(win):
    for f in win.frames:
        try:
            if f.locator("#genes_txt").count():
                return f
        except Exception:
            pass
    return None


def wait_plot(f, oid, timeout=120):
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            if f.locator(f"#{oid} img").first.get_attribute("src", timeout=500):
                return round(time.time() - t0, 1)
        except Exception:
            pass
        time.sleep(0.5)
    raise SystemExit(f"plot #{oid} did not render within {timeout}s")


exe = find_exe(unpacked)
print("launching", exe, flush=True)
proc = subprocess.Popen([exe, "--no-sandbox", "--remote-debugging-port=9333"],
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
try:
    with sync_playwright() as p:
        browser = None
        for _ in range(90):
            try:
                browser = p.chromium.connect_over_cdp("http://127.0.0.1:9333")
                break
            except Exception:
                time.sleep(1)
        if browser is None:
            raise SystemExit("could not attach to the app")
        t0 = time.time()
        f = None
        while f is None and time.time() - t0 < 300:
            for ctx in browser.contexts:
                for win in ctx.pages:
                    f = app_frame(win) or f
                    if f:
                        page = win
            time.sleep(1)
        if f is None:
            raise SystemExit("the Shiny app did not start within 300 s")
        page.set_viewport_size({"width": 1280, "height": 860})
        print(f"overview after {time.time() - t0 + wait_plot(f, 'dot', 300):.0f}s", flush=True)
        page.screenshot(path=f"{out}/1_overview.png")
        f.evaluate("document.getElementById('atlas').selectize.setValue('braun2023')")
        time.sleep(3)
        wait_plot(f, "dot")
        for tab, oid, name in [("Au fil du développement", "time", "2_time"),
                               ("Coexpression", "coex", "3_coex")]:
            f.get_by_role("tab", name=tab).click()
            print(tab, wait_plot(f, oid), "s", flush=True)
            time.sleep(1)
            page.screenshot(path=f"{out}/{name}.png")
        f.get_by_role("tab", name="Guide").click()
        time.sleep(2)
        page.screenshot(path=f"{out}/4_guide.png")
        errors = [e for e in f.locator(".shiny-output-error:not(.shiny-output-error-validation)").all_inner_texts()
                  if e.strip()]
        print("errors:", errors)
        if errors:
            raise SystemExit(1)
finally:
    proc.terminate()
