# hdid -nomount ram://2048000  # 256000 (sector) * 512 (bytes/sector)
# newfs_hfs /dev/disk6
# mkdir /tmp/ramdisk
# mount -t hfs /dev/disk6 /tmp/ramdisk
#
# umount /tmp/ramdisk/
# hdiutil detach /dev/disk6
export TMPDIR=/tmp/ramdisk
python -W ignore src/alpsdocking.py files/benzamidine.pdb files/trypsin.pdb files/trypsin_cavs.pdb files/benzamidine.itp 2 ~/best.pdb
# open ~/best.pdb
