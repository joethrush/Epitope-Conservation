import difflib
from itertools import permutations
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


ResidueKey = Tuple[str, str]


def residue_sequence(residue) -> str:
    from Bio.SeqUtils import seq1

    try:
        return seq1(residue.resname)
    except Exception:
        return "X"


def chain_sequence(chain) -> str:
    out = []
    for residue in chain:
        if residue.id[0] == " ":
            out.append(residue_sequence(residue))
    return "".join(out)


def needleman_wunsch(seq1: str, seq2: str, gap: int = -2, mismatch: int = -1,
                     match: int = 2) -> Tuple[int, List[Tuple[int, Optional[int]]]]:
    """Return a simple global alignment score and one representative."""
    n, m = len(seq1), len(seq2)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * gap
        trace[i][0] = "up"
    for j in range(1, m + 1):
        score[0][j] = j * gap
        trace[0][j] = "left"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diagonal, up, left)
            score[i][j] = best
            trace[i][j] = "diag" if diagonal == best else ("up" if up == best else "left")

    i, j = n, m
    alignment = []
    while i or j:
        direction = trace[i][j]
        if direction == "diag":
            i -= 1
            j -= 1
            alignment.append((i, j))
        elif direction == "up":
            i -= 1
            alignment.append((i, None))
        else:
            j -= 1
            alignment.append((None, j))
    alignment.reverse()

    residue_map = [(s1, s2) for s1, s2 in alignment
                   if s1 is not None and s2 is not None
                   and seq1[s1] != "-" and seq2[s2] != "-"]
    return score[n][m], residue_map


def protein_residues(chain) -> List:
    return [r for r in chain if r.id[0] == " "]


def align_chains(chain1, chain2) -> Tuple[float, Dict[int, int]]:
    residues1 = protein_residues(chain1)
    residues2 = protein_residues(chain2)
    seq1 = "".join(residue_sequence(r) for r in residues1)
    seq2 = "".join(residue_sequence(r) for r in residues2)
    score, alignment = needleman_wunsch(seq1, seq2)
    mapping = {}
    for i, j in alignment:
        mapping[i] = j
    identity = sum(residue_sequence(residues1[i]) == residue_sequence(residues2[j])
                   for i, j in mapping.items())
    coverage = identity / max(1, min(len(residues1), len(residues2)))
    return coverage, mapping


def map_chains(chains1: Sequence, chains2: Sequence) -> Tuple[Dict[str, Dict[str, str]],
                                                               List[str]]:
    """Map source chain IDs to target chain IDs, reporting ambiguous alternatives."""
    if len(chains1) == len(chains2) and len(chains1) > 0:
        mappings = list(permutations(chains2))
    else:
        mappings = [chains2]

    candidates = []
    for perm in mappings:
        mean_scores = []
        chain_map = {}
        for c1, c2 in zip(chains1, perm):
            coverage, residue_mapping = align_chains(c1, c2)
            mean_scores.append(coverage)
            chain_map[c1.id] = c2.id
        candidates.append((sum(mean_scores) / len(mean_scores), chain_map))

    best_score = max(score for score, _ in candidates)
    tied = [chain_map for score, chain_map in candidates
            if abs(score - best_score) < 1e-9]
    best = tied[0]
    alternatives = []
    if len(tied) > 1:
        alternatives.append(
            "Ambiguous chain mapping; tied alternatives: " +
            "; ".join(",".join(sorted(m.items())) for m in tied)
        )
    return best, alternatives
