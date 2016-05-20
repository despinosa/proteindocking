set TEMP=R:/TEMP
python -W ignore src/alpsdocking.py files/3ptb/benzamidine.pdb files/3ptb/trypsin.pdb files/3ptb/cavs_trypsin.pdb files/3ptb/benzamidine.itp 1 .
::python -W ignore src/alpsdocking.py files/transp_coc/cocaine.pdb files/transp_coc/transp.pdb files/transp_coc/transp_cavs.pdb files/transp_coc/cocaine.itp 1 .
::python -W ignore src/alpsdocking.py files/zinc_11909586.pdb files/1ubq.pdb files/cavities.pdb files/zinc_11909586.itp 0 .
::python -W ignore src/alpsdocking.py files/1ec3/inhibitor.pdb files/1ec3/hivhydro.pdb files/1ec3/cavs_hivhydro.pdb files/1ec3/inhibitor.itp 1 .