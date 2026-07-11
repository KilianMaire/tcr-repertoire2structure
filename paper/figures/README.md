# Figure manifest

Each figure, its message, its generating script, its data, and its status. All
figures regenerate from committed data (the structure render also needs the PyMOL
micromamba env; the two ray-traced PNGs are committed so the compositor runs
without it).

| figure | file | message | script | data | status |
| --- | --- | --- | --- | --- | --- |
| Fig 1 | `fig1_structure.png` | what a TCR-pMHC complex is (reader aid) | `plot_structure_3d.py` + `render_structure_pymol.py` | one cognate CIF in `runs/` | final |
| Fig 2 | `fig2_validation.png` | honest annotation and the leakage guard | `plot_validation.py` | `data/validation_*.csv` | final |
| Fig 3 | `fig3_retrieval.png` | structure vs sequence, discovery vs pre-registered held-out | `plot_retrieval.py` | `data/tcr_retrieval_top1.csv` | final |
| Fig 4 | `fig4_confidence_variance.png` | confidence separates TCRs, not epitopes | `plot_confidence_variance.py` | raw folds, reconstructed-only | final |
| Fig 5 | `fig5_mhc_presentation.png` | confidence reads MHC-peptide presentation | `plot_mhc_scramble.py` | raw folds, reconstructed-only | final |
| Fig 6 | `fig6_two_axis_map.png` | synthesis: presentation yes, recognition no | `plot_synthesis.py` | anchored to Fig 3/5 numbers | final |

Intermediate `_struct_view1.png` / `_struct_view2.png` are the raw PyMOL renders
that Fig 1 composes; kept so Fig 1 rebuilds without PyMOL installed.

## Significance reporting note

Retrieval significance is reported as an exact binomial test of Top-1 against naive
per-panel chance (0.25 discovery 4-way, 0.5 held-out binary), not the
label-permutation p or the TCR-blind null in the raw benchmark reports, which an
audit showed are miscalibrated on the discovery panel (they assign p near 1e-4 even
to the negative control). By the correct binomial test the discovery readouts beat
chance and the groove control does not; the held-out primary (11/18) does not
(p about 0.24), consistent with the failed pre-registration. All TCR-interface
figures use reconstructed TCRs only (poly-G stub folds excluded).
