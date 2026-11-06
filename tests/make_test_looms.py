"""Generate small synthetic .loom files that mimic the layout of the
mousebrain.org files (dev_all.loom, l5_all.loom) for testing the pipeline.
The values are random: they only exercise the code, they mean nothing."""
import sys, numpy as np, loompy

rng = np.random.default_rng(1)
out = sys.argv[1] if len(sys.argv) > 1 else "tests/fixtures"
import os; os.makedirs(out, exist_ok=True)

genes = np.array(["Pax6", "Sox2", "Eomes", "Neurod2", "Tbr1", "Gad2", "Olig2",
                  "Aqp4", "Mki67", "Hes5"] + [f"Gene{i}" for i in range(190)])

def make(path, n, ages, classes, tissues, age_numeric=False, total_attr=None):
    G = len(genes)
    cls = rng.choice(classes, n)
    m = rng.poisson(0.3, size=(G, n)).astype(np.float32)
    # give each class a marker among the first 10 genes so plots show structure
    for k, c in enumerate(classes):
        m[k % 10, cls == c] += rng.poisson(4, (cls == c).sum())
    age = rng.choice(ages, n)
    col = {
        "CellID": np.array([f"c{i}" for i in range(n)]),
        "Age": age.astype(float) if age_numeric else age.astype(str),
        "Class": cls,
        "Subclass": np.array([f"{c}_{rng.integers(1,3)}" for c in cls]),
        "Tissue": rng.choice(tissues, n),
    }
    if total_attr:
        col[total_attr] = m.sum(0)
    loompy.create(path, m, {"Gene": genes}, col)
    print("wrote", path, m.shape)

make(f"{out}/dev_test.loom", 3000, ["e11.5", "e12.5", "e13.5", "e14.0"],
     ["Radial glia", "Neuroblast", "Neuron", "Glioblast", "Vascular"],
     ["Forebrain", "Midbrain", "Hindbrain"])
make(f"{out}/dev_test_numeric.loom", 1500, ["12.5", "13.5"],
     ["Radial glia", "Neuron"], ["Forebrain"], age_numeric=True)
make(f"{out}/adol_test.loom", 3000, ["p19", "p20", "p21", "p23", "p25", "6w"],
     ["Neurons", "Astrocytes", "Oligos", "Microglia", "Vascular"],
     ["Ctx1", "Hipp", "CB", "DRG", "ENS", "SpinalC"], total_attr="TotalMolecules")
