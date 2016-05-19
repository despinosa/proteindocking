set TEMP=R:/TEMP
::python -W ignore src/alpsdocking.py files/tryp_benz/benzamidine.pdb files/tryp_benz/trypsin.pdb files/tryp_benz/trypsin_cavs.pdb files/tryp_benz/benzamidine.itp 1 .
:: python -W ignore src/alpsdocking.py files/transp_coc/cocaine.pdb files/transp_coc/transp.pdb files/transp_coc/transp_cavs.pdb files/transp_coc/cocaine.itp 1 .
::python -W ignore src/alpsdocking.py files/zinc_11909586.pdb files/1ubq.pdb files/cavities.pdb files/zinc_11909586.itp 0 .
python -W ignore src/alpsdocking.py files/1ec3/inhibitor.pdb files/1ec3/hivhydro.pdb files/1ec3/cavs_hivhydro.pdb files/1ec3/inhibitor.itp 1 .