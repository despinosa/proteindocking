hdid -nomount ram://256000  # 256000 (sector) * 512 (bytes/sector)
newfs_hfs /dev/disk1
mkdir /tmp/space
mount -t hfs /dev/disk1 /tmp/space
export TMPDIR=/tmp/ramdisk
python -W ignore src/alpsdocking.py files/benzamidine.pdb files/trypsin.pdb files/trypsin_cavs.pdb files/benzamidine.itp 2 ~/best.pdb
open ~/best.pdb
