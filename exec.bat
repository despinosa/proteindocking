set TEMP=E:/TEMP
::python -W ignore src/alpsdocking.py files/cocaine.pdb files/transp.pdb files/transp_cavs.pdb files/cocaine.itp 2 best.pdb
::python -W ignore src/alpsdocking.py files/zinc_11909586.pdb files/1ubq.pdb files/cavities.pdb files/zinc_11909586.itp 1 best.pdb
python -W ignore src/alpsdocking.py files/benzamidine.pdb files/trypsin.pdb files/trypsin_cavs.pdb files/benzamidine.itp 2 best.pdb
