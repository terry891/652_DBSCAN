start=$(date +%s)
for ((i=1;i<=5;i++))
	do
		python benchmark.py
	done
end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"
