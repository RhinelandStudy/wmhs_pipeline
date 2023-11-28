#!/bin/bash

for file in `cat /opt/fsl_binaries.txt`;
 do 
  if [ -d $file ];then
   mkdir -p /opt/$file
  else
   d=`dirname $file`
   f=`basename $file`
   mkdir -p /opt/$d
   cp $d/$f /opt/$d
  fi
done


