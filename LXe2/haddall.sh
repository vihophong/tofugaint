#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "usage: ./haddall.sh outputfiles number_of_threads"
    echo "example: ./haddall.sh LXe_tall.root 24"
else
	#hadd -f LXe_tall.root LXe_t0.root LXe_t1.root LXe_t2.root LXe_t3.root LXe_t4.root LXe_t5.root LXe_t6.root LXe_t7.root
	haddcmd="hadd -f $1 "
	for i in $(seq 0 $(expr $2 - 1))
	do
	    haddcmd+="LXe_t${i}.root "
	done
	$haddcmd
fi

