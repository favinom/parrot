rm -rf uniform*.xda
rm -rf adapted*.xda
rm -rf refined*.xdr
rm -rf first.txt
rm -rf second.txt
rm -rf third.txt
rm -rf *_in.e
rm -rf *out
rm -rf *err
rm -rf *csv
rm -rf *line*csv

if [ $1 -eq 1 ]
then
	rm -rf DiffusionOut*
	rm -rf AdvectionOut*
fi

