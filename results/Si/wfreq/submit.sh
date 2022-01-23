for i in {100..1000..100}
do
    folder='std-'$i
    cd $folder
    echo $folder
    file='std-'$i'.qsub'
    sbatch $file
    cd ..
done

# link for kin-0
for i in {100..1000..100}
do
    folder='kin-0/kin-0-'$i
    cd $folder
    echo $folder
    file='kin-0-'$i'.qsub'
    sbatch $file
    cd ../..
done

### link for mixing
#for i in {20..100..20}
#do
#    folder='kin-'$i
#    cd $folder
#    for j in {100..900..100}
#    do
#        folder='kin-'$i'-'$j
#        cd $folder 
#        echo $folder
#        file='kin-'$i'-'$j'.qsub'
#        sbatch $file
#        cd ..
#    done
#    cd ..
#done
