for n in 20 40 60 80 100 
do
    for m in {100..900..100}
    do
        there='../pw/si.save'
        here='si-'$n'-'$m'.save'
        ln -s $there $here
    done
done
