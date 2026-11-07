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


def wait_plot(f, oid, old="", timeout=90):
    """Waits until output `oid` shows an image different from `old`, then lets it settle."""
    t0 = time.time()
    while time.time() - t0 < timeout:
        s = plot_src(f, oid)
        if s and s != old:
            time.sleep(1.5)
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
        pg.evaluate(CAPTION_JS, text)
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

    say("Cette app montre <b>dans quelles cellules du cerveau</b>, et <b>à quel stade</b>, vos gènes sont exprimés.", 5)
    shot("accueil")
    say("<b>Étape 1</b> : choisir un atlas. Ici, le cerveau de souris embryonnaire, E12.5 et E13.5.", 5)
    f.locator("#atlas_card").scroll_into_view_if_needed()
    say("<b>Étape 2</b> : entrer ses gènes. Un exemple : une liste de gènes associés à l'épilepsie.", 3)
    old = plot_src(f, "dot")
    f.get_by_role("button", name="Gènes d'épilepsie").click()
    time.sleep(1.5)
    say("Puis cliquer sur <b>Afficher</b>.", 1.5)
    f.get_by_role("button", name="Afficher").click()
    wait_plot(f, "dot", old)
    shot("epilepsie_souris")
    say("Une ligne par type de cellules, une colonne par gène. <b>Gros point</b> = beaucoup de cellules expriment le gène.", 6)
    say("<b>Point clair (jaune)</b> = le gène y est exprimé fort. Ex. <i>Stxbp1</i> et <i>Kcnq2</i> : surtout dans les neurones, moins dans les neuroblastes.", 7)

    say("Changeons d'atlas : le cerveau <b>humain</b> du premier trimestre (5 à 14 semaines).", 3)
    old = plot_src(f, "dot")
    choose(f, "atlas", "braun2023")
    wait_plot(f, "dot", old)
    shot("epilepsie_humain")
    say("Les mêmes gènes, un panneau par stade. Faites défiler pour voir tous les stades.", 4)
    pg.mouse.move(900, 500)
    for _ in range(6):
        pg.mouse.wheel(0, 350)
        time.sleep(0.6)
    time.sleep(1.5)
    for _ in range(6):
        pg.mouse.wheel(0, -350)
        time.sleep(0.3)

    say("Onglet <b>Au fil du développement</b> : un gène à la fois, son évolution stade par stade.", 2)
    old = plot_src(f, "time")
    f.get_by_role("tab", name="Au fil du développement").click()
    wait_plot(f, "time", old)
    shot("evolution")
    say("Les 6 types de cellules où le gène est le plus présent sont en couleur, les autres en gris.", 6)
    old = plot_src(f, "time")
    choose(f, "time_gene", "STXBP1")
    wait_plot(f, "time", old)
    say("<i>STXBP1</i> dans les neurones : la valeur varie beaucoup d'un stade à l'autre.", 4)
    say("Chaque stade vient de peu d'embryons et de régions différentes : lisez la tendance générale, pas un point isolé.", 6)
    shot("evolution_stxbp1")

    say("Onglet <b>Coexpression</b> : deux gènes sont-ils exprimés <b>dans les mêmes cellules</b> ?", 2)
    old = plot_src(f, "coex")
    f.get_by_role("tab", name="Coexpression").click()
    wait_plot(f, "coex", old)
    old = plot_src(f, "coex")
    choose(f, "cx_stage", "11,5 sem.")
    wait_plot(f, "coex", old)
    old = plot_src(f, "coex")
    choose(f, "cx_ct", "Neuron")
    wait_plot(f, "coex", old)
    say("Ici : les neurones humains à 11,5 semaines.", 3)
    f.locator("#coex").scroll_into_view_if_needed()
    pg.mouse.move(900, 500)
    pg.mouse.wheel(0, 250)
    time.sleep(1.5)
    shot("coexpression")
    say("<b>Rouge</b> = ensemble plus souvent que le hasard (×2 = deux fois plus). <b>Bleu</b> = moins souvent. Étoiles = écart solide.", 8)
    say("Choisissez le stade et le type de cellules en haut. Le détail des paires est dans le tableau en dessous.", 5)

    pg.mouse.wheel(0, -2000)
    time.sleep(1)
    say("Onglet <b>Tableau</b> : tous les chiffres, téléchargeables pour Excel.", 2)
    f.get_by_role("tab", name="Tableau").click()
    time.sleep(3)
    shot("tableau")
    time.sleep(3)

    say("Onglet <b>Guide</b> : les quatre atlas, ce que veulent dire les chiffres, et les précautions.", 2)
    f.get_by_role("tab", name="Guide").click()
    time.sleep(2)
    shot("guide")
    time.sleep(4)
    say("À vous : collez votre liste de gènes et cliquez sur <b>Afficher</b>.", 4)
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
