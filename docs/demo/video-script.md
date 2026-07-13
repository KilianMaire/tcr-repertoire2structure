# R2S submission video: narration and shot list

Target length 3 minutes (hard max). Narration is about 420 words, roughly 145 words per
minute. Researcher track. The visuals are captured by `scripts/capture_demo.py`, which replays
the real run_1 transcript into the web app so every frame is genuine agent output.

## Narration (what you say)

**[0:00 to 0:20] The problem.**
A researcher sequences a T cell repertoire and gets tens of thousands of TCRs. They want to
know two things: what do these T cells recognise, and what do the recognising complexes look
like. Today that means stitching together a specificity database, a germline reconstructor, a
folding model, and a lot of caution. Nobody hands the researcher one honest answer.

**[0:20 to 0:45] What R2S is.**
Repertoire2Structure is that missing tool. You drop a 10x contig CSV in the browser, and a team
of Claude agents takes it from raw repertoire to QC'd TCR to peptide MHC structures, with honest
specificity annotation and a skeptical check on every structure. It was built entirely with
Claude Code and the Claude Agent SDK.

**[0:45 to 1:15] Intake and annotation.**
An intake agent interviews you: the data, your question, where to compute. Then it hands off.
Annotation is by similarity, never prediction. Each clonotype is matched by TCRdist against
labelled references, and it is annotated only when a neighbour is close enough, always with the
distance and a confidence tier. Here, of a hundred and nineteen clonotypes, twenty three are
annotated with high confidence and ninety are flagged unannotatable. No label is ever forced.

**[1:15 to 1:55] Routing and folding.**
A structure strategist reasons about the question and routes each group of clonotypes to the
right tool. This is a class one TCR to peptide MHC question, so it stays on Protenix, the default
workhorse; other questions route to TCRdock or MHC Fine. An executor builds the fold artifact.
And crucially, before a single fold is spent, any clonotype whose V domain could not be
reconstructed is flagged as a stub, so you are never sold a structure built on a placeholder.

**[1:55 to 2:35] The honesty that is the point.**
A predicted structure does not confirm specificity. Folding models impose canonical docking
geometry even on non binding sequences, so the QC agent is a skeptic: it calibrates each fold
against a scrambled peptide control and refuses to call a structure trustworthy when it cannot
beat that control. Our preregistered benchmark found the same thing a parallel study found this
month: structural confidence reads peptide presentation, not TCR recognition. We do not oversell
it.

**[2:35 to 3:00] The contribution.**
So the contribution is not a new model. It is an agentic, honest pipeline from a biological
question to a structure, with the guardrails a researcher actually needs: leakage aware
annotation, stub flagging, and scramble calibrated QC. It is open source. Drop a repertoire, and
Claude gives you an honest map from sequence to structure, without touching a pipeline.

## Shot list (what is on screen, matched to the narration)

| Time | On screen | Captured file |
|------|-----------|---------------|
| 0:00 | Hero landing page, dropzone with the DNA icon | 01_landing.png |
| 0:20 | CSV dropped, run starts, first agent line appears | 02_run_start.png |
| 0:45 | Intake agent question, then the green "your turn" indicator and the chat box | 03_intake_your_turn.png |
| 1:00 | annotate_specificity result: 23 high, 90 unannotatable | 04_annotation.png |
| 1:15 | structure strategist routing the group to Protenix | 05_routing.png |
| 1:35 | executor building the fold notebook (and the stub warning if any) | 06_build_artifact.png |
| 2:00 | full timeline scrolled, tool pills and agent lanes visible | 07_timeline_full.png |
| 2:35 | sidebar with the finished run, "agent working" to idle | 08_done.png |

The full screen recording is saved as `demo_capture/session.webm`; use it as b roll under the
narration, or re record the same flow live with your own voice.
