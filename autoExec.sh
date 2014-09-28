for countOuter in  1 2 3 4 5 6 7 8 9 10
do
	for  count in 1 2 3 4 5 6 7 8 9 10
	do
		echo " $countOuter . $count -> executing... "
		sh both.sh >> firstSetup.txt
		echo "------------------" >> firstSetup.txt

	done
done