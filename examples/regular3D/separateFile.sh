filename=$1;
lne=$2;

intro=17;
ne=`expr ${lne} \\* ${lne} \\* ${lne}`;
nn=`expr ${lne} + 1`;
nn=`expr ${nn} \\* ${nn} \\* ${nn}`;

rm -rf first.txt
rm -rf second.txt
rm -rf third.txt

line=`expr ${intro} + ${ne} + ${nn} + 1`

tail -n +${line} ${filename} > third.txt

line=`expr ${intro} + ${ne} + ${nn}`

head -n ${line} ${filename} > temp.txt

line=`expr ${intro} + ${ne} + 1`

tail -n +${line} temp.txt > second.txt

line=`expr ${intro} + ${ne}`

head -n ${line} temp.txt > first.txt

rm -rf temp.txt
