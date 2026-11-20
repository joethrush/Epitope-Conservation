# Epitope Conservation Tool - Design Summary

## Goal

Build a Python command-line tool that compares two structural states of a
protein and ranks conserved, accessible epitope candidates. The first use case
is pre-fusion versus post-fusion viral glycoproteins, but the workflow should be
usable for other proteins.

## Inputs and chain selection

- Accept one PDB or mmCIF file for each state.
- Let the user designate one or more target chains in each structure. Candidate
  epitopes may span those chains.
- Retain other chains as structural context for accessibility and occlusion
  calculations.
- Use biological assemblies when they can be determined from file metadata.
  If the assembly cannot be determined, warn and require an explicit override
  before analyzing the provided chains as-is.
- Match chains and residues by sequence alignment, not by PDB chain labels or
  residue numbering. Evaluate plausible ambiguous chain mappings and report
  alternatives rather than silently choosing one.

## Structure preparation

- Leave the input files unchanged.
- Write cleaned structures and analysis outputs to a separate output directory.
- Remove unrelated heteroatoms while preserving modeled glycans for occlusion
  analysis.
- Identify N-linked glycosylation sequons from the protein sequence and report
  modeled glycans found in the structures. Do not predict O-linked
  glycosylation.
- For alternate atom locations, use the conformer with the highest occupancy,
  with a deterministic tie-break.
- Analyze every model in a multi-model structure.

## Candidate epitope search

1. In each structure, consider each residue in the designated target chains as
   a seed.
2. Form an eight-residue neighborhood from the seed and its seven nearest
   eligible C-alpha atoms within 8 angstroms. Neighbors may be on other
   designated target chains.
3. Omit a seed if it has fewer than seven eligible neighbors within that
   distance.
4. Compare the same eight sequence-mapped residues in the other structure.
   All eight residues must map and have the required backbone coordinates in
   both structures. They do not have to remain nearest neighbors in the other
   structure.
5. Merge identical residue-set neighborhoods and equivalent symmetry-related
   chain copies. Keep overlapping but distinct neighborhoods as separate
   candidates.

## Structural comparison and accessibility

- Fit each candidate independently using its mapped backbone atoms (N, C-alpha,
  C, and O) and calculate its backbone RMSD. A global whole-protein fit is not
  used to determine local epitope geometry.
- Calculate neighborhood preservation as the number of mapped residues shared
  by the two eight-residue neighborhoods, including the seed, divided by 8.
- Evaluate accessibility in each full assembly, including other chains and
  modeled glycans. Accessibility is a structural proxy, not a guarantee of
  antibody binding.
- A candidate must have at least five of its eight residues with relative
  solvent accessibility of at least 20% in each structure.
- Report N-linked sequons and modeled glycan locations separately, including
  chain and residue identifiers.

## Ranking

Return the top 100 candidates. Use these provisional, documented score
components:

| Component | Weight | Initial calculation |
| --- | ---: | --- |
| Backbone geometry | 60% | `max(0, 1 - RMSD / 4 angstroms)` |
| Neighborhood preservation | 20% | Shared mapped neighborhood residues divided by 8 |
| Surface exposure | 10% | The lower of the two structures' mean relative solvent accessibility values |
| Glycan occlusion | 10% | The lower fraction of patch surface remaining exposed after modeled glycans are included |

Show component scores alongside the combined score. If glycan data are
incomplete for either structure, mark the glycan component unavailable and
renormalize the other weights rather than assuming no occlusion.

For multiple models, compare every model from the first input with every model
from the second. Rank candidates using the median composite score across model
pairs and report the score range. Break exact ranking ties by lower backbone
RMSD, then higher exposure, then residue and chain identifiers.

## Outputs

- A ranked CSV containing candidate residues and chain IDs, combined and
  component scores, RMSD, exposure in each state, glycan impact, and model
  variability.
- A separate glycosylation report for N-linked sequons and modeled glycans.
- Cleaned structures and a PyMOL script for visualizing candidate patches.
- A separate mapping report when ambiguous chain assignments affect results.

## Initial scope notes

- The initial candidate size is eight residues; support for changing this size
  can be added later.
- The 8-angstrom neighborhood distance, 4-angstrom RMSD scale, and score
  weights are initial defaults. Score weights should be configurable and may
  later be calibrated against known examples.
- Accessibility is assessed structurally; antibody-scale approach geometry is
  outside the initial scope.
