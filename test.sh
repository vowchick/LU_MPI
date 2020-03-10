for ((k = 1; k <= 4; k++)) do
  for ((i = 0; i <= 30; i++)) do
    for ((j = 3; j <= $i; j+=3)) do
mpirun -np $k ./a.out $i $j 
done;
done;
done;
