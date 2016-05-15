export TMPDIR=/dev/shm/pifi
python -W ignore src/alpsdocking.py files/benzamidine.pdb files/trypsin.pdb files/trypsin_cavs.pdb files/benzamidine.itp 1 ~
# open ~/best.pdb
