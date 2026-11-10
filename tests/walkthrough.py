"""Walks through the app like a first-time user, with French captions, and records
a video plus one screenshot per step. Also serves as the end-to-end browser test:
it fails if any Shiny output shows an error.

Usage: python tests/walkthrough.py <app_url> <out_dir>
Writes <out_dir>/walkthrough.webm and <out_dir>/NN_<step>.png.
"""
import os, shutil, sys, time
from playwright.sync_api import sync_playwright

url, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
W, H = 1440, 900

CAPTION_JS = """(t) => {
  let c = document.getElementById('wt-caption');
  if (!c) {
    c = document.createElement('div'); c.id = 'wt-caption';
    c.style.cssText = 'position:fixed;left:50%;bottom:28px;transform:translateX(-50%);z-index:99999;' +
      'max-width:1100px;background:rgba(20,40,60,.92);color:#fff;font:500 22px/1.4 system-ui,sans-serif;' +
      'padding:14px 24px;border-radius:10px;box-shadow:0 4px 18px rgba(0,0,0,.3);text-align:center;';
    document.body.appendChild(c);
  }
  c.style.display = t ? 'block' : 'none'; c.innerHTML = t;
}"""


def app_frame(pg):
    for f in pg.frames:
        try:
            if f.locator("#genes_txt").count():
                return f
        except Exception:
            pass
    return None


def plot_src(f, oid):
    try:
        return f.locator(f"#{oid} img").first.get_attribute("src", timeout=500) or ""
    except Exception:
        return ""


def choose(f, input_id, value):
    """Selects a value in a Shiny selectInput through its visible selectize widget."""
    f.evaluate("([id, v]) => document.getElementById(id).selectize.setValue(v)", [input_id, value])


def busy(f, oid):
    """True while Shiny recomputes output `oid` or shows a progress bar."""
    try:
        return f.evaluate(
            "(id) => { const e = document.getElementById(id);"
            " return !!(e && e.classList.contains('recalculating')) ||"
            " document.querySelectorAll('.shiny-notification').length > 0; }", oid)
    except Exception:
        return True


def wait_plot(f, oid, old="", timeout=90):
    """Waits until output `oid` shows a new image and Shiny has finished computing."""
    t0 = time.time()
    while time.time() - t0 < timeout:
        s = plot_src(f, oid)
        if s and s != old and not busy(f, oid):
            time.sleep(1.5)
            if not busy(f, oid):
                return time.time() - t0
        time.sleep(0.3)
    raise RuntimeError(f"plot #{oid} did not render within {timeout}s")


with sync_playwright() as p:
    b = p.chromium.launch()
    ctx = b.new_context(viewport={"width": W, "height": H},
                        record_video_dir=out, record_video_size={"width": W, "height": H})
    pg = ctx.new_page()
    step = [0]

    def say(text, hold):
        # On a first visit shinylive reloads the page once its service worker is
        # installed; a caption drawn during that reload is retried on the new page.
        for _ in range(20):
            try:
                pg.evaluate(CAPTION_JS, text)
                break
            except Exception:
                time.sleep(0.5)
        time.sleep(hold)

    def shot(name):
        step[0] += 1
        pg.screenshot(path=f"{out}/{step[0]:02d}_{name}.png")

    pg.goto(url)
    say("Chargement de l'app… (une dizaine de secondes la première fois)", 0.5)
    f = None
    t0 = time.time()
    while f is None and time.time() - t0 < 240:
        f = app_frame(pg)
        time.sleep(0.5)
    if f is None:
        raise RuntimeError("app did not load")
    wait_plot(f, "dot", timeout=240)
    print(f"app ready after {time.time() - t0:.0f}s", flush=True)

    say("Cette app montre <b>dans quelles cellules du cerveau</b>, <b>dans quelle région</b> et <b>à quel stade</b> vos gènes sont exprimés.", 5)
    shot("accueil")
    say("<b>Étape 1</b> : l'atlas. Celui-ci réunit 8 études : le prosencéphale de souris de E10.5 à l'adulte, région par région.", 6)
    say("Les âges sont regroupés autour de E10.5, E12.5, E13.5, E15.5, P0, P10, P21, P40 et P100. Le signe ≈ indique l'âge réel le plus proche.", 6)
    say("<b>Étape 2</b> : les gènes. Ici <i>Myo9b, Trio, Rhoa, Rac1, Cdc42</i>, déjà affichés : on clique sur <b>Afficher</b> après chaque changement.", 4)
    shot("vue_ensemble")
    say("Une ligne par type cellulaire, une colonne par gène, un panneau par âge. <b>Gros point</b> = beaucoup de cellules expriment le gène ; <b>jaune</b> = exprimé fort.", 7)
    say("Les types cellulaires sont les mêmes pour toutes les études : progéniteurs Nkx2-1+, interneurones Lhx6+, PV+, SST+, pyramidaux, glie…", 6)
    pg.mouse.move(900, 500)
    for _ in range(5):
        pg.mouse.wheel(0, 400)
        time.sleep(0.6)
    time.sleep(1)
    for _ in range(5):
        pg.mouse.wheel(0, -400)
        time.sleep(0.3)

    say("<b>Étape 3</b> : affiner. Choisir des <b>régions</b> : par exemple la MGE et le cortex somatosensoriel.", 3)
    old = plot_src(f, "dot")
    choose(f, "regions", ["MGE", "Cortex somatosensoriel"])
    wait_plot(f, "dot", old)
    shot("regions")
    say("Seules les cellules de ces régions restent. La MGE n'existe qu'à l'embryon : les panneaux postnataux montrent le cortex.", 6)

    say("Et ne garder que les cellules qui expriment un gène, ici <i>Lhx6</i> : la lignée des interneurones issus de la MGE.", 3)
    f.locator("#pos_gene").fill("Lhx6")
    old = plot_src(f, "dot")
    f.get_by_role("button", name="Afficher").click()
    wait_plot(f, "dot", old)
    shot("lhx6")
    say("Vos gènes, uniquement dans les cellules <i>Lhx6</i>+, de la MGE embryonnaire jusqu'au cortex.", 6)

    say("Onglet <b>Au fil du développement</b> : un gène à la fois, âge par âge.", 2)
    old = plot_src(f, "time")
    f.get_by_role("tab", name="Au fil du développement").click()
    wait_plot(f, "time", old)
    old = plot_src(f, "time")
    choose(f, "time_gene", "Trio")
    wait_plot(f, "time", old)
    shot("evolution_trio")
    say("<i>Trio</i> dans les cellules <i>Lhx6</i>+ : chaque ligne est un type cellulaire. Un point manque quand le groupe a moins de 20 cellules.", 7)
    say("Les âges viennent d'études différentes : lisez la tendance générale, pas un point isolé.", 5)

    f.locator("#pos_gene").fill("")
    choose(f, "regions", [])
    old = plot_src(f, "coex")
    f.get_by_role("button", name="Afficher").click()
    say("Onglet <b>Coexpression</b> : deux gènes sont-ils exprimés <b>dans les mêmes cellules</b> ?", 2)
    f.get_by_role("tab", name="Coexpression").click()
    wait_plot(f, "coex", old)
    old = plot_src(f, "coex")
    choose(f, "cx_stage", "≈ E13.5")
    wait_plot(f, "coex", old)
    old = plot_src(f, "coex")
    choose(f, "cx_ct", "Interneurones Lhx6+ immatures (issus de la MGE)")
    wait_plot(f, "coex", old)
    f.locator("#coex").scroll_into_view_if_needed()
    pg.mouse.move(900, 500)
    pg.mouse.wheel(0, 250)
    time.sleep(1.5)
    shot("coexpression")
    say("Ici : les interneurones Lhx6+ vers E13.5. <b>Rouge</b> = ensemble plus souvent que le hasard, <b>bleu</b> = moins souvent, étoiles = écart solide.", 8)
    pg.mouse.wheel(0, -2000)
    time.sleep(1)

    say("Toute une famille de gènes : <b>Ou charger une liste prête</b> → les 159 régulateurs des Rho GTPases.", 3)
    choose(f, "preset", "rho:all")
    time.sleep(1.5)
    f.get_by_role("tab", name="Carte de chaleur").click()
    old = plot_src(f, "heat")
    f.get_by_role("button", name="Afficher").click()
    wait_plot(f, "heat", old, timeout=180)
    shot("carte_rho")
    say("Onglet <b>Carte de chaleur</b> : une ligne par gène, rangés par famille (Rho GTPases, GEF, GAP, GDI, myosines IX). <b>Foncé</b> = beaucoup de cellules.", 7)
    pg.mouse.move(900, 500)
    for _ in range(6):
        pg.mouse.wheel(0, 450)
        time.sleep(0.6)
    time.sleep(1)
    pg.mouse.wheel(0, -4000)
    time.sleep(1)

    say("Onglet <b>Tableau</b> : tous les chiffres, à télécharger pour Excel.", 2)
    f.get_by_role("tab", name="Tableau").click()
    time.sleep(3)
    shot("tableau")
    say("Onglet <b>Guide</b> : les études, les âges réels, ce que veulent dire les chiffres et les précautions.", 2)
    f.get_by_role("tab", name="Guide").click()
    time.sleep(2)
    shot("guide")
    time.sleep(4)
    say("À vous : choisissez un atlas, collez vos gènes et cliquez sur <b>Afficher</b>.", 4)
    say("", 0.5)

    errors = []
    for fr in pg.frames:
        try:
            errors += [e for e in fr.locator(".shiny-output-error:not(.shiny-output-error-validation)").all_inner_texts() if e.strip()]
        except Exception:
            pass
    video = pg.video.path()
    ctx.close()
    b.close()
    shutil.move(video, f"{out}/walkthrough.webm")
    print("errors:", errors)
    sys.exit(1 if errors else 0)
