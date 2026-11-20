from pathlib import Path
from typing import Iterable, List, Sequence
import numpy as np

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


GLYCAN_NAMES = {
    "NAG", "NDG", "MAN", "BMA", "GAL", "SIA", "FUC", "GLC", "GCU",
    "XYP", "XYL", "A2G", "GlcNAc", "Man", "Gal",
}


def load_structure(path: str):
    """Parse a PDB or mmCIF file and return its Bio.PDB structure object."""
    path = str(path)
    if path.lower().endswith(".cif") or path.lower().endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
    return parser.get_structure(Path(path).stem, path)


def ca_coordinate(residue):
    """Return the residue's selected C-alpha coordinates, or None if absent."""
    if "CA" not in residue:
        return None
    atom = highest_occupancy_altloc(residue["CA"])
    return atom.coord if atom is not None else None


def highest_occupancy_altloc(atom):
    """Return the atom object used for coordinates without further altloc selection."""
    if atom.altloc in (" ", "A"):
        return atom
    return atom


def clean_structure(structure):
    """Remove unrecognized hetero residues and select occupied atom conformers."""
    for model in structure:
        for chain in list(model):
            for residue in list(chain):
                if residue.id[0] == " ":
                    continue
                hetname = residue.resname.strip().upper()
                if residue.id[0] != " " and hetname not in GLYCAN_NAMES:
                    chain.detach_child(residue.id)
                else:
                    for atom in residue:
                        conformers = [atom]
                        if atom.is_disordered():
                            conformers = list(atom.disordered_get_list())
                        if atom.is_disordered():
                            chosen = max(
                                conformers,
                                key=lambda a: (a.occupancy if a.occupancy is not None else 0.0,
                                                -ord(a.altloc if a.altloc not in (" ", "") else "A"))
                            )
                            atom.disordered_select(chosen.altloc)
    return structure


def write_cleaned(structure, path: str):
    """Write a Bio.PDB structure to a PDB file at the requested path."""
    io = PDBIO()
    io.set_structure(structure)
    io.save(path)


def get_chain(structure, model_id: int, chain_id: str):
    """Return a named chain from a model, or None when the chain is absent."""
    model = structure[model_id]
    return model[chain_id] if chain_id in model else None


def assembly_metadata_available(path: str) -> bool:
    """Check PDB or mmCIF metadata for biological assembly annotations."""
    path = str(path)
    if path.lower().endswith((".cif", ".mmcif")):
        try:
            data = MMCIF2Dict(path)
            return bool(data.get("_pdbx_struct_assembly_gen"))
        except Exception:
            return False
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return any(line.startswith("REMARK 350") for line in handle)
    except OSError:
        return False
