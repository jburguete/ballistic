#!/bin/bash
sed -i "s/"$1"\."$2"\."$3"/"$4"\."$5"\."$6"/g" $1.$2.$3/Doxyfile README.md
git mv $1.$2.$3 $4.$5.$6
