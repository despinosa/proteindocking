# hdid -nomount ram://2048000  # 256000 (sector) * 512 (bytes/sector)
# newfs_hfs /dev/disk6
# mkdir /tmp/ramdisk
# mount -t hfs /dev/disk6 /tmp/ramdisk
#
# umount /tmp/ramdisk/
# hdiutil detach /dev/disk6
# export TMPDIR=/tmp/ramdisk
python -W ignore src/alpsdocking.py files/tryp_benz/benzamidine.pdb files/tryp_benz/trypsin.pdb files/tryp_benz/trypsin_cavs.pdb files/tryp_benz/benzamidine.itp 1 .
# open $TMPDIR/best.pdb
