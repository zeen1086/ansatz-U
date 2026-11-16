"""Quality diagnostics must not become production admission gates."""
from contextlib import redirect_stderr
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'7_spin'))
import common
from complete_readout import metric


class QualityDiagnosticsTests(unittest.TestCase):
    def test_regression_is_optional_even_when_failed_or_invalid(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(common, 'SPINS', Path(tmp)):
            report = Path(tmp)/'regression/RESULT.json'
            report.parent.mkdir()
            for payload in (None, '{', '[]', json.dumps({'passed': False}),
                            json.dumps({'passed': True, 'model_id': 'different'})):
                if payload is not None:
                    report.write_text(payload)
                log = io.StringIO()
                with redirect_stderr(log):
                    common.check_regression()
                self.assertIn('QUALITY WARNING', log.getvalue())
            report.write_text(json.dumps({'passed': True, 'model_id': common.CONFIG['model_id']}))
            with redirect_stderr(io.StringIO()) as log:
                common.check_regression()
            self.assertEqual(log.getvalue(), '')

    def test_small_positive_overlap_runs_but_singular_overlap_fails(self):
        with redirect_stderr(io.StringIO()) as log:
            s, t, _ = metric(np.eye(30)[None]*.1)
        self.assertIn('QUALITY WARNING', log.getvalue())
        np.testing.assert_allclose(t@s@t, np.eye(60), atol=1e-12)
        with self.assertRaisesRegex(ValueError, 'Nonpositive'):
            metric(np.zeros((1,30,30)))

    def test_nonfinite_spectrum_still_fails(self):
        with self.assertRaisesRegex(ValueError, 'Invalid spectra'):
            common.chemical_potential(np.array([[np.nan, 1.]]), np.ones(1), 1)


if __name__ == '__main__':
    unittest.main()
