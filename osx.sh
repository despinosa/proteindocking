if ["$(mount | grep /tmp/ramdisk)" = ""]; then
    export RAMDEV=$(hdid -nomount ram://2048000)  # 256000 (sector) * 512 (bytes/sector)
    newfs_hfs $RAMDEV
    mkdir /tmp/ramdisk
    mount -t hfs $RAMDEV /tmp/ramdisk
fi
export TMPDIR=/tmp/ramdisk
python -W ignore src/alpsdocking.py files/tryp_benz/benzamidine.pdb files/tryp_benz/trypsin.pdb files/tryp_benz/trypsin_cavs.pdb files/tryp_benz/benzamidine.itp 1 .
open ./best*.pdb &
umount /tmp/ramdisk/
hdiutil detach $RAMDEV
