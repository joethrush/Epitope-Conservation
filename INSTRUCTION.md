# Epitope Conservation Pipeline

## Purpose

Compare two PDB or mmCIF protein structures and rank conserved, accessible
eight-residue epitope patches. It can be used to compare pre-fusion and
post-fusion states or other pairs of related protein structures.

## Installation

Run from the project directory:

```powershell
python -m pip install -e .
```

## Required inputs

- A first PDB or mmCIF structure.
- A second PDB or mmCIF structure.
- Target chain IDs in the first structure.
- Target chain IDs in the second structure.

Chain IDs are supplied separately because chain names or residue numbering may
differ between structures.

## Basic operation

```powershell
python -m epitope_conservation.cli state1.pdb state2.pdb `
  --chains1 A,B `
  --chains2 A,B `
  --force-as-is `
  --output results
```

Use `--force-as-is` when the chains present in the supplied files are the
intended unit of analysis. The implementation can detect PDB/mmCIF assembly
metadata, but does not yet expand biological assemblies.

## Options

| Option | Default | Meaning |
|---|---:|---|
| `--chains1` | Required | Comma-separated target chains in the first structure |
| `--chains2` | Required | Comma-separated target chains in the second structure |
| `--output`, `-o` | `epitope_output` | Output directory |
| `--neighborhood-size` | `8` | Residues in each candidate patch |
| `--distance` | `8.0` | C-alpha neighbor cutoff in angstroms |
| `--rmsd-scale` | `4.0` | RMSD scale used for geometry score |
| `--weight-backbone` | `0.60` | Backbone geometry weight |
| `--weight-neighborhood` | `0.20` | Neighborhood preservation weight |
| `--weight-exposure` | `0.10` | Solvent exposure weight |
| `--weight-glycan` | `0.10` | Glycan occlusion weight |
| `--force-as-is` | Off | Analyze the supplied chains without assembly expansion |

View the same options with:

```powershell
python -m epitope_conservation.cli --help
```

## Analysis workflow

1. Parse both PDB/mmCIF structures.
2. Remove unrelated heteroatoms while retaining modeled glycan residues.
3. Write cleaned PDB copies to the output directory.
4. Match chains and residues by sequence alignment.
5. Build patches from each eligible target-chain residue and its nearest C-alpha
   neighbors within the configured distance.
6. Map each patch to the corresponding residues in the second structure.
7. Require all required backbone atoms: N, C-alpha, C, and O.
8. Independently fit each mapped patch and calculate local backbone RMSD.
9. Estimate solvent accessibility and reject patches with fewer than five
   residues above 20% accessibility in either structure.
10. Apply geometry, neighborhood, exposure, and glycan scores.
11. Compare every model pair for multi-model structures.
12. Rank candidates using median composite score across model pairs and return
    the top 100.

## Scoring

The default combined score is:

```text
0.60 * backbone geometry
+ 0.20 * neighborhood preservation
+ 0.10 * lower surface exposure
+ 0.10 * lower glycan-adjusted exposure
```

- **Backbone geometry:** `max(0, 1 - RMSD / 4)` by default.
- **Neighborhood preservation:** fraction of mapped patch residues preserved in
  the second-state neighborhood.
- **Surface exposure:** lower mean relative exposure of the two states.
- **Glycan occlusion:** lower fraction of patch residues remaining distant from
  modeled glycan atoms.

If glycan scoring is unavailable, the remaining component weights are
renormalized.

Candidate accessibility is a structural proxy, not evidence that an antibody
will bind.

## Outputs

All outputs are written beneath `--output`.

- `ranked_candidates.csv`: primary top-100 result table.
- `glycosylation_report.csv`: predicted N-linked sequons and modeled glycans.
- `mapping_report.txt`: chain mappings used for each model pair.
- `candidate_patches.pml`: PyMOL script for visualizing candidate patches.
- `state1_cleaned.pdb`: cleaned first structure.
- `state2_cleaned.pdb`: cleaned second structure.
- `warnings.txt`: assembly/chain-as-is warning, when applicable.

### Ranked candidate columns

`ranked_candidates.csv` includes:

- Rank and residue/chain identifiers.
- Median and range of composite score.
- Backbone, neighborhood, exposure, and glycan component scores.
- Median RMSD and RMSD range.
- Median exposure and per-residue RSA values for both states.
- Number of model pairs used.
- Glycan availability status.

### Visualize candidates in PyMOL

```text
@candidate_patches.pml
```

Run this command from the directory containing the cleaned PDB files and PML
script.

## Current limitations

- Biological assemblies are detected but not expanded.
- Chain mapping selects the best pairing rather than reporting every plausible
  alternative.
- Modeled-glycan and accessibility calculations are structural approximations.
- Alternate-atom handling is only partially implemented.
- N-linked sequons are predicted from protein sequence; O-linked glycosylation
  is not predicted.
