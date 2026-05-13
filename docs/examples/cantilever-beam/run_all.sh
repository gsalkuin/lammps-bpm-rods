#!/bin/bash

# File to store PIDs
PID_FILE="run_pids.txt"
> "$PID_FILE"   # clear PID file at start

# === LOOP OVER FILES ===
for f in lam/*.lam; do
    # extract number between N and .lam
    fname=$(basename "$f")
    num=$(echo "$fname" | sed -n 's/.*-\([0-9]\+\)\.lam/\1/p')
    
    target="$num"

    # make target folder
    rm -rf $target
    mkdir $target

    # copy files
    cp in.bend-point $target/
    cp $f $target/

    # run command in background, with log
    (
        cd $target || exit
        sed -i "s/NBONDS/$num/g" in.bend-point
        nohup ~/Data/lammps-dev-20250710/lmp_bpm-sym -in in.bend-point > output.log 2>&1 &
        echo $! >> ../"$PID_FILE"
    )
done

wait   # wait for all background jobs to finish

