if ["$(mount | grep /tmp/ramdisk)" = ""]; then
    export RAMDEV=$(hdid -nomount ram://3072000)  # 256000 (sector) * 512 (bytes/sector)
    newfs_hfs $RAMDEV
    mkdir /tmp/ramdisk
    mount -t hfs $RAMDEV /tmp/ramdisk
fi
export TMPDIR=/tmp/ramdisk
python -W ignore src/alpsdocking.py files/3ptb/benzamidine.pdb files/3ptb/trypsin.pdb files/3ptb/cavs_trypsin.pdb files/3ptb/benzamidine.itp 1 .
open ./best*.pdb &
# umount /tmp/ramdisk/
# hdiutil detach $RAMDEV
