# Epitope Conservation

Rank conserved, accessible epitope candidates between two PDB or mmCIF protein structures.

## Install

```powershell
python -m pip install -e .
```

## Run

```powershell
python -m epitope_conservation.cli state1.pdb state2.pdb `
  --chains1 A,B --chains2 A,B --force-as-is `
  --output epitope_output
```

Use `--force-as-is` only when you intentionally want to analyze the chains in the
provided files without biological-assembly expansion.

The pipeline writes:

- `ranked_candidates.csv`
- `glycosylation_report.csv`
- `mapping_report.txt`
- `candidate_patches.pml`
- `state1_cleaned.pdb`
- `state2_cleaned.pdb`

Default neighborhood, RMSD, and ranking settings are documented in
`DESIGN_SUMMARY.md`. Override them with `--distance`, `--neighborhood-size`,
`--rmsd-scale`, and the four `--weight-*` options.
