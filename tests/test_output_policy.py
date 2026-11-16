"""Compact outputs must retain restart inputs until every native result is complete."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import json
import subprocess
import unittest
from unittest.mock import patch

import local_vfunc


class OutputPolicyTests(unittest.TestCase):
    def test_failed_or_truncated_hamiltonian_keeps_inputs(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);model='a'*64;basis='b'*64
            background=root/'pure.au';background.write_bytes(b'background')
            Path(str(background)+'.json').write_text(json.dumps(
                dict(kind='pure_H0_no_local_fields_no_HU',background_id='fixture-background')))
            rvr=root/'B_general'
            rvr.write_bytes(b'LANIO2_LRVR_V2  '+bytes(32)+model.encode()+basis.encode())
            (root/'PROJECTION.json').write_text(json.dumps(dict(basis_id=basis,background_id='fixture-background')))

            def evaluate(directory,structure,out,sites,step,**options):
                out.mkdir()
                (out/'value.bin').write_bytes(b'coefficients')
                (out/'dsite_0001.bin').write_bytes(b'derivatives')
                (out/'PREDICTION.json').write_text('{}')
                return dict(status='COMPLETE',model_id=model,L=1)

            for failure in ('process','truncated'):
                with self.subTest(failure=failure):
                    out=root/failure
                    args=SimpleNamespace(model=root,structure=root/'structure.npz',rvr=rvr,
                        background=background,kpoints=root/'k.dat',binary=root/'cal_Hsk',
                        output=out,derivative_sites=[1],step=.001,components=False,keep_intermediates=False)

                    def native(command,**options):
                        label=Path(command[1]).stem
                        if label!='value' and failure=='process':
                            raise subprocess.CalledProcessError(1,['cal_Hsk'])
                        (out/f'{label}.index').write_text('\n'.join([
                            'LANIO2_LOCAL_H_V2','1 60 1 F',model,basis,'1 0 0 0 1 0 0 1'])+'\n')
                        (out/f'{label}.k000001.h').write_bytes(bytes(16*60**2) if label=='value' else b'truncated')

                    with patch.object(local_vfunc,'metadata',return_value={'model_id':model}), \
                         patch.object(local_vfunc,'evaluate',side_effect=evaluate), \
                         patch.object(local_vfunc.subprocess,'run',side_effect=native):
                        with self.assertRaises((subprocess.CalledProcessError,RuntimeError)):
                            local_vfunc.hamiltonian(args)
                    self.assertTrue((out/'value.bin').is_file())
                    self.assertTrue((out/'value.nml').is_file())
                    self.assertTrue((out/'dsite_0001.bin').is_file())
                    self.assertTrue((out/'value.index').is_file())
                    self.assertTrue((out/'PREDICTION.json').is_file())
                    self.assertFalse((out/'HAMILTONIAN.json').exists())


if __name__=='__main__':unittest.main()
