from __future__ import annotations

def annotation_metrics(annotations, labels) -> dict:
    n = len(annotations)
    annotated = [a for a in annotations if a.annotatable]
    correct = sum(1 for a in annotated if labels.get(a.clonotype_id) == a.epitope)
    n_labeled = sum(1 for a in annotations if a.clonotype_id in labels)
    unann = sum(1 for a in annotations if not a.annotatable)
    precision = correct / len(annotated) if annotated else float("nan")
    recall = correct / n_labeled if n_labeled else float("nan")
    return {
        "precision": precision,
        "recall": recall,
        "unannotatable_rate": unann / n if n else float("nan"),
        "n": n,
        "n_annotated": len(annotated),
        "n_correct": correct,
    }
