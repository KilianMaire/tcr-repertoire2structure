# Figure manifest

Each figure, its message, its generating script, and its data. Every figure uses
the shared style module `scripts/figstyle.py` (Okabe-Ito palette, no titles, no
result text on the plot, bold panel tags, PDF plus PNG output) and regenerates from
committed data with no GPU. Structure panels embed committed PyMOL renders, so the
compositors run without the PyMOL env.

## Graphical abstract

`graphical_abstract` (`plot_graphical_abstract.py`): the one self-contained visual.
Input structure, then the two specificity questions, presentation answered (cognate
vs scramble groove by per-residue pLDDT, AUROC up to 0.99) and recognition not
(held-out 0.61, 11/18, p 0.24), then the one-line takeaway. Being caption-free, it
is the only figure that carries the headline finding and key numbers as text.

## Main figures (dense multi-panel)

| figure | file | message | script | data |
| --- | --- | --- | --- | --- |
| Fig 1 | `fig1_structure` | pipeline and an annotated TCR-pMHC complex (5 panels: pipeline, whole complex, groove, interface, two axes) | `plot_structure_3d.py` + `render_structure_pymol.py` | committed renders `_struct_view1/2.png`, `_interface_flu.png` |
| Fig 2 | `fig2_validation` | annotation and the leakage guard (4 panels) | `plot_validation.py` | `data/validation_*.csv`, `docs/validation_donor1_metrics.json` |
| Fig 3 | `fig3_retrieval` | structure vs sequence: discovery, held-out, controls, binomial null (5 panels) | `plot_retrieval.py` | `data/tcr_retrieval_top1.csv` |
| Fig 4 | `fig4_confidence_variance` | confidence separates TCRs, not epitopes (4 panels) | `plot_confidence_variance.py` | `data/confidence_variance.csv` |
| Fig 5 | `fig5_mhc_presentation` | confidence reads MHC-peptide presentation (5 panels, incl. groove structure) | `plot_mhc_scramble.py` | `data/mhc_presentation.csv`, `data/peptide_plddt.csv`, `data/scramble_anchor_permissiveness.csv`, committed renders |
| Fig 6 | `fig6_two_axis_map` | synthesis: presentation yes, recognition no | `plot_synthesis.py` | anchored to Fig 3/5 numbers |

## Supplementary figures

| figure | file | message | script | data |
| --- | --- | --- | --- | --- |
| Fig S1 | `figS1_retrieval_strata` | full retrieval battery, all folds vs reconstructed | `plot_supp_retrieval_strata.py` | `data/tcr_retrieval_top1.csv` |
| Fig S2 | `figS2_groove_confidence` | per-residue confidence, cognate vs scramble | `plot_supp_groove_conf.py` | committed renders `_groove_conf_*.png` |
| Fig S3 | `figS3_chain_pair_iptm` | chain-pair interface confidence matrix | `plot_supp_chain_pair.py` | `data/chain_pair_iptm_example.csv` |
| Fig S4 | `figS4_reproducibility` | per-sample reproducibility of the readouts | `plot_supp_reproducibility.py` | `data/per_sample_readouts.csv` |
| Fig S5 | `figS5_complex_gallery` | gallery of complexes across epitopes, one canonical TCR-up orientation | `plot_supp_gallery.py` (composits) + `render_gallery.py` (selects and renders) | committed renders `_gallery{1..4}_*.png` |

Committed intermediate renders (raw PyMOL PNGs prefixed `_`) let the structure
figures rebuild without PyMOL installed. To regenerate a render, run
`scripts/render_structure_pymol.py` in the pymol micromamba env (modes: complex,
groove, groove_conf, interface, single, single_uniform). The gallery renders are
produced by `scripts/render_gallery.py`, which per epitope selects the best
reconstructed cognate clonotype (highest TCR-to-peptide ipTM, poly-G stubs excluded)
from the run manifests and renders it in the shared TCR-up frame, so the panel
provenance is recorded rather than ad hoc. Per-sample and chain-pair CSVs come from
`scripts/extract_supp_data.py` (run once, needs `runs/`).

## Significance reporting note

Retrieval significance is an exact binomial test of Top-1 against naive per-panel
chance (0.25 discovery, 0.5 held-out), not the label-permutation p or TCR-blind null
in the raw benchmark reports, which an audit showed are miscalibrated on the
discovery panel (they assign p near 1e-4 even to the negative control). By the
binomial test the discovery readouts beat chance and the groove control does not;
the held-out primary (11/18) does not (p about 0.24), consistent with the failed
pre-registration. All TCR-interface figures use reconstructed TCRs only (poly-G stub
folds excluded).
