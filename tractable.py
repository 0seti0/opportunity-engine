"""Covalent-tractability gate — drop ANTIBODY-CLASS / biologic-only targets (cell-surface antigens,
secreted proteins, integrins, circulating cytokines) whose druggable site is extracellular, so a covalent
small molecule can't engage a productive site.

This is a CLASS/localization gate, NOT a cysteine gate: a ligandable-cysteine test gives false negatives
(NLRP3 Cys279, PARP7, STAT6 are real covalent targets whose cysteines are scout-fragment-invisible) and
false positives (MS4A1/CD20 has a stray Cys but is a pure antibody target). Intracellular enzymes,
kinases, transporters, nuclear receptors and TFs all PASS; only the curated biologic classes are dropped.
"""

BIOLOGIC = {
    # immune checkpoints + T/B-cell surface antigens (antibody / bispecific / CAR-T / ADC)
    "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VSIR", "BTLA", "ICOS", "CD27", "CD28",
    "TNFRSF4", "TNFRSF9", "TNFRSF18", "CD40", "CD40LG", "CD80", "CD86", "SLAMF7", "KLRG1", "CD96",
    "CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD5", "CD7", "CD8A", "CD19", "MS4A1", "CD22", "CD23", "FCER2",
    "CD24", "CD33", "CD37", "CD38", "CD44", "CD47", "SIRPA", "CD52", "CD70", "CD74", "CD79A", "CD79B",
    "CD200", "SDC1", "IL2RA", "IL3RA", "NCAM1", "CD276", "VTCN1",
    # tumor surface antigens
    "EPCAM", "TACSTD2", "MSLN", "FOLR1", "FOLH1", "CEACAM5", "MUC1", "MUC16", "DLL3", "NECTIN4", "CLDN18",
    "GPC3", "TNFRSF17", "GUCY2C", "NPR1", "NPR2",
    # integrins
    "ITGA4", "ITGAL", "ITGAV", "ITGB2", "ITGB3", "ITGB6", "ITGB7", "ITGA2B",
    # secreted / coagulation / complement / circulating
    "TNF", "LTA", "IL1B", "IL2", "IL4", "IL5", "IL6", "IL13", "IL17A", "IL23A", "IL31", "TSLP", "IGHE",
    "VEGFA", "ANGPT2", "PCSK9", "ANGPTL3", "APOC3", "LPA", "SOST", "DKK1", "INHBE", "GDF15",
    "F2", "F5", "F7", "F9", "F10", "F11", "F12", "F13A1", "FGA", "FGB", "FGG", "VWF", "PROC", "KLKB1",
    "PLG", "ALB", "SERPINC1", "SERPIND1", "SERPINA1", "C3", "C5", "C1QA", "MASP2", "CFB", "CFD",
    # cytokine/GF receptors hit EXTRACELLULARLY by biologics (NOT the RTKs, which keep their kinase domain)
    "IL4R", "IL5RA", "IL6R", "IL23R", "IL31RA", "CRLF2", "TNFRSF1A", "TNFRSF1B", "FCGRT",
}


def covalently_tractable(sym):
    """True unless the target is a curated antibody-class / secreted / biologic-only target."""
    return sym not in BIOLOGIC


if __name__ == "__main__":
    for s in ("CD274", "CD19", "MS4A1", "NLRP3", "TIPARP", "STAT6", "RAF1", "ERBB2", "FGFR3", "SLC1A5"):
        print(f"  {s:8} {'PASS (covalently tractable)' if covalently_tractable(s) else 'DROP (antibody-class/biologic)'}")
