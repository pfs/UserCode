#! /bin/sh

base=$(basename $1)
echo "base $base"
if [ -d $2 ];
then
    echo "Begin file copy. File 1 $1 -> file 2 $2/$base"
else
    echo "Begin file copy. File 1 $1 -> file 2 $2"
fi

ls -al $1

if [ ! -f $1 ];
then
    echo "NO FILE TO COPY"
    exit 1
fi
# if [ -f $2 ];
# then
#     eos root://eosuser.cern.ch// rm $2
# fi

xrdcp -f -C adler32 $1 root://eosuser.cern.ch/$2
res=$?
checksum1=`md5sum $1 | awk '{print $1;}'`
checksum2=0
if [ -d $2 ]; 
then
    checksum2=`md5sum $2/${base} | awk '{print $1;}'`
else
    checksum2=`md5sum $2 | awk '{print $1;}'`
fi
echo "checksum1 $checksum1 checksum2 $checksum2"
while [ $res -ne 0 ];
#"$checksum1" != "$checksum2" ];
do
    echo "Files differ"
    # if [ -f $2 ];
    # then
    # 	eos root://eosuser.cern.ch// rm $2
    # fi
    xrdcp -f -C adler32 $1 root://eosuser.cern.ch/$2
    res=$?
    echo "Copy repeated"
    checksum1=`md5sum $1 | awk '{print $1;}'`
    checksum2=0
    if [ -d root://eosuser.cern.ch/$2 ]; 
    then
	base=$(basename $1)
	cheksum2=`md5sum $2/${base} | awk '{print $1;}'`
    else
	cheksum2=`md5sum $2 | awk '{print $1;}'`
    fi
      
    echo "checksum1 $checksum1 checksum2 $checksum2"

done
echo "File copy done"
if [ -d $2 ]; 
then
    base=$(basename $1)
    eos root://eosuser.cern.ch// ls -al $2/$base
else
    :
#    eos root://eosuser.cern.ch// ls -al $2/$base
fi

rm $1
