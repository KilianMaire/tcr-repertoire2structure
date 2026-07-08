# Fold procedure (real fold_fn)

The Claude agent implements fold_fn by driving Protenix on Google Colab through Playwright, following the established procedure: open the Colab notebook, upload the construct FASTA and precomputed MSA, run `protenix_base_default_v1.0.0` at 5 seeds, let background execution survive disconnect, and download the resulting CIF models to `~/.playwright-mcp/`. It returns the local model paths. The loop is resumable: a job with a `fold_<id>.done.txt` marker is skipped.
