for n in {20..100..20}
do
  for m in {100..900..100}
  do
    qsub='si-'$n'-'$m'.qsub'
    echo $qsub
    sbatch $qsub
  done
done
