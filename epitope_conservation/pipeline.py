from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import csv
import statistics
import numpy as np
from scipy.spatial.distance import cdist

from .mapping import align_chains, chain_sequence, protein_residues, residue_sequence
from .structure import (load_structure, clean_structure, write_cleaned, ca_coordinate,
                        GLYCAN_NAMES, assembly_metadata_available)


BACKBONE_ATOMS = ("N", "CA", "C", "O")


@dataclass(frozen=True)
class ResidueRef:
    """Store a residue's chain, structure identifier, and three-letter name."""

    chain_id: str
    residue_id: tuple
    residue_name: str

    @property
    def label(self):
        """Format the residue reference for reports and candidate labels."""
        return f"{self.chain_id}:{self.residue_name}{self.residue_id[1]}{self.residue_id[2].strip()}"


@dataclass
class Candidate:
    """Represent an epitope patch by its residue labels and structure keys."""

    residues: List[ResidueRef]
    mapped_residue_keys: List[Tuple[str, tuple]]


def residue_ref(residue):
    """Convert a Bio.PDB residue into the compact reference used by the pipeline."""
    return ResidueRef(residue.parent.id, residue.id, residue.resname)


def backbone_coords(residue):
    """Return N, CA, C, and O coordinates, or None if any backbone atom is absent."""
    coords = []
    for name in BACKBONE_ATOMS:
        if name not in residue:
            return None
        atom = residue[name]
        coords.append(atom.coord)
    return np.asarray(coords, dtype=float)


def kabsch_rmsd(coords1, coords2):
    """Superimpose two coordinate sets and return their root-mean-square deviation."""
    coords1 = np.asarray(coords1, dtype=float)
    coords2 = np.asarray(coords2, dtype=float)
    c1 = coords1.mean(axis=0)
    c2 = coords2.mean(axis=0)
    x = coords1 - c1
    y = coords2 - c2
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    fitted = y @ rotation
    return float(np.sqrt(np.mean(np.sum((x - fitted) ** 2, axis=1))))


def residue_atoms(structure, selected_chain_ids):
    """Collect protein atoms, glycan atoms, and their residues from chosen chains."""
    atoms = []
    glycans = []
    residues = []
    glycan_residues = []
    for model in structure:
        for chain in model:
            if chain.id not in selected_chain_ids:
                continue
            for residue in chain:
                if residue.id[0] == " ":
                    residues.append(residue)
                    atoms.extend([atom for atom in residue.get_atoms() if atom.element != "H"])
                elif residue.resname.strip().upper() in GLYCAN_NAMES:
                    glycans.extend([atom for atom in residue.get_atoms() if atom.element != "H"])
                    glycan_residues.append(residue)
    return atoms, glycans, residues, glycan_residues


def point_surface_accessibility(residue, all_atoms, probe_radius=1.4, n_points=96):
    """Estimate residue exposure by sampling points around each heavy atom.

    Points farther than the atom radius from every structure atom are counted as
    exposed; the returned value is the exposed fraction of sampled points.
    """
    ca = ca_coordinate(residue)
    if ca is None:
        return 0.0
    coords = np.asarray([atom.coord for atom in all_atoms], dtype=float)
    residue_atom_coords = np.asarray([atom.coord for atom in residue.get_atoms()
                                      if atom.element != "H"], dtype=float)
    if residue_atom_coords.size == 0:
        return 0.0

    z = np.linspace(-1.0, 1.0, 12)
    phi = np.pi * (3.0 - np.sqrt(5.0))
    points = []
    for i, zi in enumerate(z):
        radius = np.sqrt(max(0.0, 1.0 - zi * zi))
        theta = i * phi
        points.append((radius * np.cos(theta), radius * np.sin(theta), zi))
    points = np.asarray(points, dtype=float)

    exposed = 0
    for atom in residue_atom_coords:
        sphere = atom + (probe_radius + 1.7) * points
        distances = cdist(sphere, coords)
        exposed += int(np.count_nonzero(np.min(distances, axis=1) > 1.7))
    return exposed / max(1, len(residue_atom_coords) * len(points))


def residue_rsa(residue, all_atoms):
    """Scale the estimated exposure to a relative accessibility score from 0 to 1."""
    # The compact 1.7-A probe proxy is scaled to a documented relative scale.
    accessibility = point_surface_accessibility(residue, all_atoms)
    max_accessibility = 0.60
    return min(1.0, accessibility / max_accessibility)


def make_candidates(model, target_chains: Sequence[str], chain_map: Optional[Dict[str, str]] = None,
                    size: int = 8, distance: float = 8.0):
    """Build unique fixed-size residue patches around C-alpha neighbor groups.

    Each protein residue with a C-alpha atom acts as a seed. Its closest eligible
    residues within ``distance`` form a candidate, and duplicate patches are
    removed before returning.
    """
    chains = [model[c] for c in target_chains if c in model]
    residues = [r for chain in chains for r in protein_residues(chain)]
    ca_coords = {}
    for residue in residues:
        coord = ca_coordinate(residue)
        if coord is not None:
            ca_coords[(residue.parent.id, residue.id)] = np.asarray(coord)
    candidates = []
    for seed in residues:
        seed_key = (seed.parent.id, seed.id)
        if seed_key not in ca_coords:
            continue
        neighbors = []
        for residue in residues:
            key = (residue.parent.id, residue.id)
            if key not in ca_coords:
                continue
            if key == seed_key:
                continue
            d = float(np.linalg.norm(ca_coords[seed_key] - ca_coords[key]))
            if d <= distance:
                neighbors.append((d, key, residue))
        neighbors.sort(key=lambda x: (x[0], x[1][0], str(x[1][1])))
        if len(neighbors) < size - 1:
            continue
        neighborhood = [(0.0, seed_key, seed)] + neighbors[:size - 1]
        refs = [residue_ref(r) for _, _, r in neighborhood]
        keys = [(r.parent.id, r.id) for _, _, r in neighborhood]
        candidates.append(Candidate(refs, keys))

    unique = {}
    for candidate in candidates:
        key = tuple(sorted((c.chain_id, c.residue_id) for c in candidate.residues))
        unique.setdefault(key, candidate)
    return list(unique.values())


def map_candidate(candidate, residue_map):
    """Map every candidate residue to the other structure, failing if one is unmapped."""
    mapped = []
    for residue in candidate.residues:
        key = (residue.chain_id, residue.residue_id)
        target = residue_map.get(key)
        if target is None:
            return None
        mapped.append(target)
    return mapped


def build_residue_mappings(model1, model2, chains1, chains2):
    """Map each available source chain to its best-aligned target chain and residues."""
    mapping = {}
    chain_pairs = {}
    for c1 in chains1:
        if c1 not in model1:
            continue
        best = None
        best_coverage = -1.0
        for c2 in chains2:
            if c2 not in model2:
                continue
            coverage, residue_mapping = align_chains(model1[c1], model2[c2])
            if coverage > best_coverage:
                best_coverage = coverage
                best = (c2, residue_mapping)
        if best is not None:
            c2, residue_mapping = best
            chain_pairs[c1] = c2
            r1 = protein_residues(model1[c1])
            r2 = protein_residues(model2[c2])
            for i, j in residue_mapping:
                if i < len(r1) and j < len(r2):
                    mapping[(c1, r1[i].id)] = (c2, r2[j].id)
    return mapping, chain_pairs


def glycan_component(candidate_mapped_residues, model1, model2, all_glycans1, all_glycans2,
                     structure_ready=True):
    """Return the candidate's fraction of backbone points clear of source glycans.

    Return None when the structure is not ready or no candidate backbone
    coordinates are available.
    """
    if not structure_ready:
        return None
    coords1 = []
    coords2 = []
    for residue in candidate_mapped_residues:
        coords = backbone_coords(residue)
        if coords is not None:
            coords1.append(coords)
    if not coords1:
        return None
    coords1 = np.vstack(coords1)
    glycan1 = np.asarray([atom.coord for atom in all_glycans1], dtype=float)
    if len(glycan1) == 0:
        fraction1 = 1.0
    else:
        fraction1 = float(np.mean(np.min(cdist(coords1, glycan1), axis=1) > 5.0))

    mapped_source_count = len(candidate_mapped_residues)
    return fraction1


def patch_glycan_exposure(residues, glycans):
    """Measure the fraction of patch backbone coordinates over 5 A from glycans."""
    coords = []
    for residue in residues:
        coords_value = backbone_coords(residue)
        if coords_value is not None:
            coords.append(coords_value)
    if not coords:
        return None
    patch = np.vstack(coords)
    if len(glycans) == 0:
        return 1.0
    glycan_coords = np.asarray([atom.coord for atom in glycans], dtype=float)
    distances = cdist(patch, glycan_coords)
    return float(np.mean(np.min(distances, axis=1) > 5.0))


def all_model_atoms(model):
    """Split a model's non-hydrogen atoms into protein and recognized glycan atoms."""
    atoms = []
    glycans = []
    for chain in model:
        for residue in chain:
            if residue.id[0] == " ":
                atoms.extend([atom for atom in residue.get_atoms() if atom.element != "H"])
            elif residue.resname.strip().upper() in GLYCAN_NAMES:
                glycans.extend([atom for atom in residue.get_atoms() if atom.element != "H"])
    return atoms, glycans


def sequons(chain):
    """Find N-X-S/T glycosylation motifs and return their asparagine residues."""
    residues = protein_residues(chain)
    sequence = chain_sequence(chain)
    residue_map = {i: residue for i, residue in enumerate(residues)}
    found = []
    for i in range(len(sequence) - 2):
        if sequence[i] == "N" and sequence[i + 1] != "P" and sequence[i + 2] in ("S", "T"):
            residue = residue_map.get(i)
            if residue is not None:
                found.append(residue_ref(residue))
    return found


def modeled_glycans(model):
    """Return references to recognized glycan residues already present in a model."""
    found = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " " and residue.resname.strip().upper() in GLYCAN_NAMES:
                found.append(residue_ref(residue))
    return found


def run_pipeline(state1, state2, chains1, chains2, output_dir="epitope_output",
                 neighborhood_size=8, distance=8.0, rmsd_scale=4.0,
                 weights=(0.60, 0.20, 0.10, 0.10), force_as_is=False):
    """Analyze two structures, rank matching epitope patches, and write reports.

    The workflow validates assembly metadata, loads and cleans both structures,
    aligns chains and candidate patches, calculates component scores, aggregates
    results across model pairs, then writes CSV, text, and PyMOL output files.
    """
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    warnings = []

    assembly1 = assembly_metadata_available(state1)
    assembly2 = assembly_metadata_available(state2)
    if (assembly1 or assembly2) and not force_as_is:
        raise ValueError(
            "Biological assembly metadata is available but automatic assembly expansion "
            "is not enabled for this input. Pass --force-as-is only if analyzing "
            "the provided chains as-is is intentional."
        )
    if not assembly1 or not assembly2:
        warnings.append(
            "No biological assembly metadata was detected for at least one input; "
            "analysis used the provided chains as-is."
        )

    structure1 = clean_structure(load_structure(state1))
    structure2 = clean_structure(load_structure(state2))

    write_cleaned(structure1, str(output / "state1_cleaned.pdb"))
    write_cleaned(structure2, str(output / "state2_cleaned.pdb"))

    results = []
    mapping_report = []
    all_candidate_keys = set()

    # Compare every model pair so alternate models contribute to the final ranking.
    for model1 in structure1:
        if model1.id not in structure1:
            continue
        candidates1 = make_candidates(model1, chains1, size=neighborhood_size, distance=distance)
        for model2 in structure2:
            candidates2 = make_candidates(model2, chains2, size=neighborhood_size, distance=distance)
            mapping, chain_pairs = build_residue_mappings(model1, model2, chains1, chains2)
            if chain_pairs:
                mapping_report.append(
                    f"Model {model1.id}->{model2.id}: " +
                    ", ".join(f"{source}->{target}" for source, target in sorted(chain_pairs.items()))
                )

            atoms1, glycans1 = all_model_atoms(model1)
            atoms2, glycans2 = all_model_atoms(model2)

            # Keep only complete residue mappings with usable backbones and exposure.
            for candidate1 in candidates1:
                mapped_keys = map_candidate(candidate1, mapping)
                if mapped_keys is None or len(mapped_keys) != neighborhood_size:
                    continue
                source_residues = []
                target_residues = []
                valid = True
                for ref in candidate1.residues:
                    if ref.chain_id not in model1:
                        valid = False
                        break
                    source_residue = model1[ref.chain_id][ref.residue_id]
                    source_residues.append(source_residue)
                for chain_id, residue_id in mapped_keys:
                    if chain_id not in model2:
                        valid = False
                        break
                    target_residue = model2[chain_id][residue_id]
                    target_residues.append(target_residue)
                if not valid:
                    continue

                backbone1 = [backbone_coords(r) for r in source_residues]
                backbone2 = [backbone_coords(r) for r in target_residues]
                if any(v is None for v in backbone1 + backbone2):
                    continue
                rmsd = kabsch_rmsd(np.vstack(backbone1), np.vstack(backbone2))

                rsa1 = [residue_rsa(r, atoms1) for r in source_residues]
                rsa2 = [residue_rsa(r, atoms2) for r in target_residues]
                if sum(v >= 0.20 for v in rsa1) < 5 or sum(v >= 0.20 for v in rsa2) < 5:
                    continue

                source_keys = {(r.chain_id, r.residue_id) for r in candidate1.residues}
                target_key_set = set(mapped_keys)
                target_seed_key = mapped_keys[0]
                target_residue_neighbors = []
                target_seed_residue = model2[target_seed_key[0]][target_seed_key[1]]
                seed_coord = ca_coordinate(target_seed_residue)
                if seed_coord is not None:
                    for target_chain_id in chains2:
                        if target_chain_id not in model2:
                            continue
                        for residue in protein_residues(model2[target_chain_id]):
                            coord = ca_coordinate(residue)
                            if coord is None:
                                continue
                            delta = float(np.linalg.norm(seed_coord - np.asarray(coord)))
                            if delta <= distance:
                                target_residue_neighbors.append(
                                    (delta, (residue.parent.id, residue.id))
                                )
                target_residue_neighbors.sort(key=lambda item: (
                    item[0], item[1][0], str(item[1][1])
                ))
                target_neighborhood = {key for _, key in target_residue_neighbors[:neighborhood_size]}
                neighborhood_score = len(target_key_set & target_neighborhood) / neighborhood_size
                exposure = min(float(np.mean(rsa1)), float(np.mean(rsa2)))
                glycan1 = patch_glycan_exposure(source_residues, glycans1)
                glycan2 = patch_glycan_exposure(target_residues, glycans2)
                glycan_score = None if glycan1 is None or glycan2 is None else min(glycan1, glycan2)

                # Combine the measured features, renormalizing if glycan data is absent.
                backbone_score = max(0.0, 1.0 - rmsd / rmsd_scale)
                if glycan_score is None:
                    active_weights = (weights[0], weights[1], weights[2])
                    active_total = sum(active_weights)
                    combined = (
                        backbone_score * active_weights[0] +
                        neighborhood_score * active_weights[1] +
                        exposure * active_weights[2]
                    ) / active_total
                else:
                    combined = (
                        backbone_score * weights[0] +
                        neighborhood_score * weights[1] +
                        exposure * weights[2] +
                        glycan_score * weights[3]
                    ) / sum(weights)

                residue_labels = "|".join(ref.label for ref in candidate1.residues)
                result_key = tuple(sorted(source_keys))
                all_candidate_keys.add(result_key)
                results.append({
                    "key": result_key,
                    "labels": residue_labels,
                    "chain_ids": ",".join(sorted(set(ref.chain_id for ref in candidate1.residues))),
                    "composite": combined,
                    "backbone": backbone_score,
                    "neighborhood": neighborhood_score,
                    "exposure": exposure,
                    "glycan": glycan_score,
                    "rmsd": rmsd,
                    "rsa1": ";".join(f"{v:.3f}" for v in rsa1),
                    "rsa2": ";".join(f"{v:.3f}" for v in rsa2),
                    "model1": model1.id,
                    "model2": model2.id,
                })

    # Aggregate repeated observations of the same patch across model pairs.
    grouped = {}
    for result in results:
        grouped.setdefault(result["key"], []).append(result)

    ranked = []
    for key, variants in grouped.items():
        composites = [v["composite"] for v in variants]
        rmsds = [v["rmsd"] for v in variants]
        exposures = [v["exposure"] for v in variants]
        glycan_values = [v["glycan"] for v in variants if v["glycan"] is not None]
        first = variants[0]
        ranked.append({
            **first,
            "composite": statistics.median(composites),
            "composite_min": min(composites),
            "composite_max": max(composites),
            "rmsd_median": statistics.median(rmsds),
            "rmsd_min": min(rmsds),
            "rmsd_max": max(rmsds),
            "exposure_median": statistics.median(exposures),
            "glycan": statistics.median(glycan_values) if glycan_values else None,
            "model_pairs": len(variants),
        })

    # Rank by score first, then use structural fit and exposure to break ties.
    ranked.sort(key=lambda r: (
        -r["composite"],
        r["rmsd_median"],
        -r["exposure_median"],
        r["labels"],
    ))
    ranked = ranked[:100]

    # Write ranked score details and per-state glycosylation information.
    csv_path = output / "ranked_candidates.csv"
    fieldnames = [
        "rank", "residue_ids", "chain_ids", "combined_score", "backbone_score",
        "neighborhood_score", "exposure_score", "glycan_score", "rmsd_median",
        "rmsd_range", "exposure_median", "composite_range", "model_pairs",
        "rsa_state1", "rsa_state2", "glycan_status",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for rank, result in enumerate(ranked, start=1):
            glycan_status = "available" if result["glycan"] is not None else "unavailable"
            writer.writerow({
                "rank": rank,
                "residue_ids": result["labels"],
                "chain_ids": result["chain_ids"],
                "combined_score": f"{result['composite']:.6f}",
                "backbone_score": f"{result['backbone']:.6f}",
                "neighborhood_score": f"{result['neighborhood']:.6f}",
                "exposure_score": f"{result['exposure_median']:.6f}",
                "glycan_score": "" if result["glycan"] is None else f"{result['glycan']:.6f}",
                "rmsd_median": f"{result['rmsd_median']:.6f}",
                "rmsd_range": f"{result['rmsd_min']:.6f}-{result['rmsd_max']:.6f}",
                "exposure_median": f"{result['exposure_median']:.6f}",
                "composite_range": f"{result['composite_min']:.6f}-{result['composite_max']:.6f}",
                "model_pairs": result["model_pairs"],
                "rsa_state1": result["rsa1"],
                "rsa_state2": result["rsa2"],
                "glycan_status": glycan_status,
            })

    glycan_path = output / "glycosylation_report.csv"
    with open(glycan_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["state", "record_type", "chain_id", "residue_id", "residue_name"])
        for state_name, structure in (("state1", structure1), ("state2", structure2)):
            for model in structure:
                for chain in model:
                    for ref in sequons(chain):
                        writer.writerow([state_name, "N-linked sequon", ref.chain_id,
                                         str(ref.residue_id[1]), ref.residue_name])
                for ref in modeled_glycans(model):
                    writer.writerow([state_name, "modeled glycan", ref.chain_id,
                                     str(ref.residue_id[1]), ref.residue_name])

    # Save the chain mapping summary separately for easier inspection.
    mapping_path = output / "mapping_report.txt"
    with open(mapping_path, "w", encoding="utf-8") as handle:
        handle.write("Chain mapping report\n")
        handle.write("\n".join(mapping_report) if mapping_report else "No ambiguities detected.")

    # Generate a PyMOL script that highlights each ranked patch.
    pymol_path = output / "candidate_patches.pml"
    with open(pymol_path, "w", encoding="utf-8") as handle:
        handle.write("load state1_cleaned.pdb, state1\n")
        handle.write("load state2_cleaned.pdb, state2\n")
        handle.write("hide everything\n")
        handle.write("show cartoon, state1 or state2\n")
        for rank, result in enumerate(ranked, start=1):
            labels = result["labels"].split("|")
            # CSV labels are chain:idresidue, PyMOL selections use chain and residue number.
            selections = []
            for label in labels:
                chain_and_residue = label.split(":", 1)
                if len(chain_and_residue) == 2:
                    chain_id, residue_name_and_number = chain_and_residue
                    digits = "".join(c for c in residue_name_and_number if c.isdigit())
                    if digits:
                        selections.append(f"chain {chain_id} and resi {digits}")
            if selections:
                selection = " or ".join(selections)
                handle.write(f"select candidate{rank}, {selection}\n")
                handle.write(f"show sticks, candidate{rank}\n")

    if warnings:
        (output / "warnings.txt").write_text("\n".join(warnings), encoding="utf-8")

    print(f"Wrote {len(ranked)} ranked candidates to {csv_path}")
    print(f"Wrote glycosylation report to {glycan_path}")
    print(f"Wrote PyMOL script to {pymol_path}")
    return ranked
