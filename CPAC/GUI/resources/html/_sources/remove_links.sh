#!/bin/bash

echo "#!/usr/bin/env python" >> remove_links.py
echo "# -*- coding: utf-8 -*-" >> remove_links.py
echo "import fileinput"  >> remove_links.py
echo ""  >> remove_links.py

for doc in *.txt
do
    if [ $(egrep -c -o \`.*\<.*\`_ ${doc}) -eq 0 ];
    then
        continue
    fi
    echo 'with open("'$(pwd)/${doc}'", "r") as f:' >> remove_links.py
        echo -e '\tlines = f.readlines()' >> remove_links.py
        echo -e '\tnewlines = []' >> remove_links.py
        echo -e '\tfor line in lines:' >> remove_links.py
    egrep -o \`.*\<.*\`_ ${doc} | while read -r line
    do 
        echo -e '\t\tline = line.replace("'${line}'", "'$(echo ${line} | cut -d\< -f1)'".replace("`",""))' >> remove_links.py
    done
    echo -e '\t\tnewlines.append(line)' >> remove_links.py
    echo ""  >> remove_links.py
    echo 'with open("'$(pwd)/${doc}'", "w") as f:' >> remove_links.py
        echo -e '\tfor line in newlines:' >> remove_links.py
        echo -e '\t\tf.write(line)' >> remove_links.py
    echo "" >> remove_links.py
done
python remove_links.py
rm remove_links.py
