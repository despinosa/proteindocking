if ["$(mount | grep /tmp/ramdisk)" = ""]; then
    export RAMDEV=$(hdid -nomount ram://2048000)  # 256000 (sector) * 512 (bytes/sector)
    newfs_hfs $RAMDEV
    mkdir /tmp/ramdisk
    mount -t hfs $RAMDEV /tmp/ramdisk
fi
export TMPDIR=/tmp/ramdisk
python -W ignore src/alpsdocking.py files/5byy/benzodiazepin.pdb files/5byy/kinase.pdb files/5byy/cavs_kinase.pdb files/5byy/benzodiazepin.itp 1 .
open ./best*.pdb &
# umount /tmp/ramdisk/
# hdiutil detach $RAMDEV
