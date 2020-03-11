for ((k = 1; k <= 4; k++)) do
  for ((i = 12; i <= 30; i+=3)) do
    for ((j = 3; j <= 3; j+=3)) do
mpirun -np $k ./a.out $i $j 
done;
done;
done;
