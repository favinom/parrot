rm -rf uniform*.xda
rm -rf adapted*.xda
rm -rf refined*.xdr
rm -rf first.txt
rm -rf second.txt
rm -rf third.txt

if [ $1 -eq 1 ]
then
	rm -rf DiffusionOut*
	rm -rf AdvectionOut*
fi

