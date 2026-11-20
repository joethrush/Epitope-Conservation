import numpy as np
import unittest
from epitope_conservation.pipeline import kabsch_rmsd
from epitope_conservation.mapping import needleman_wunsch


class PipelineTests(unittest.TestCase):
    """Check core structural fitting and sequence-alignment behaviors."""

    def test_kabsch_identical_structures(self):
        """Identical point sets should have effectively zero aligned RMSD."""
        coordinates = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        self.assertLess(kabsch_rmsd(coordinates, coordinates), 1e-6)

    def test_kabsch_reports_rotation_difference(self):
        """Distinct point geometry should produce a positive RMSD."""
        coordinates1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        coordinates2 = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        self.assertGreater(kabsch_rmsd(coordinates1, coordinates2), 0.0)

    def test_sequence_alignment_maps_matching_residues(self):
        """An exact sequence match should align each residue at the same index."""
        score, alignment = needleman_wunsch("ACDE", "ACDE")
        self.assertGreater(score, 0)
        self.assertEqual(alignment, [(0, 0), (1, 1), (2, 2), (3, 3)])
