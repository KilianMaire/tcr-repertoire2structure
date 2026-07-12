# Table 1. Datasets

| dataset | role | unit | labels | size used | accession |
| --- | --- | --- | --- | --- | --- |
| 10x 4 donor CD8 dextramer | validation (ground truth) | paired clonotype | per-cell binarized dextramer specificity | donor 1, 3,325 labelled clonotypes | 10x Genomics public dataset |
| retrieval panels (derived) | structure vs sequence test | folded TCR-pMHC construct | dextramer cognate | A*02:01 n=48 (29 reconstructed); A*11:01 n=24 (18 reconstructed) | subset of the above |
| TABLO | application (scale) | paired clonotype | none | full repertoire, end to end | Zenodo 10.5281/zenodo.13119615 |

Notes. The retrieval panels are drawn from the validation set, restricted to
unannotatable (sequence-novel by TCRdist) clonotypes. Reconstructed counts exclude
the poly-G stub constructs (see Methods, stub contamination). The A*11:01 set is
the pre-registered held-out arm.
