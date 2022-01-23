for i in {100..1000..100}
do
    folder='std-'$i
    cd $folder
    rm *.save
    ln -s ../../pw/si.save si.save
    
    there='../../wstat/si_1000.wstat.save'
    here='si_1000.wstat.save'
    ln -s ${there} ${here}
    cd ..
done

## link for kin-0
for i in {100..1000..100}
do
    folder='kin-0/kin-0-'$i
    cd $folder
    echo $folder
    rm *.save
    ln -s ../../../pw/si.save si.save

    there='../../../wstat/si_kin_1000.wstat.save'
    here='si_kin_1000.wstat.save'
    ln -s ${there} ${here}
    cd ../..
done

## link for mixing
for i in {20..100..20}
do
    folder='kin-'$i
    cd $folder
    for j in {100..900..100}
    do
        folder='kin-'$i'-'$j
        cd $folder 
        echo $folder
        rm *.save
        ln -s ../../../pw/si.save si.save
        
        there='../../../wstat/si_kin_'$i'_'$j'.wstat.save'
        here='si_kin_'$i'_'$j'.wstat.save'
        ln -s ${there} ${here}
        cd ..
    done
    cd ..
done
