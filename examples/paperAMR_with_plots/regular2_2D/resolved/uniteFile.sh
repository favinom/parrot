filename=$1;

cat first.txt > ${filename}
cat second.txt >> ${filename}
cat third.txt >> ${filename}

rm -rf first.txt
rm -rf second.txt
rm -rf third.txt
