
i=0
for nh in 9.99999978e-03   1.41000003e-02   1.99999996e-02 \
         2.82000005e-02   3.97999994e-02   5.62000014e-02 \
         7.94000030e-02   1.12000003e-01   1.58000007e-01 \
         2.24000007e-01   3.16000015e-01   4.46999997e-01 \
         6.30999982e-01   8.90999973e-01   1.25999999e+00 \
         1.77999997e+00   2.50999999e+00   3.54999995e+00 \
         5.01000023e+00   7.07999992e+00   1.00000000e+01 \
         1.41000004e+01   2.00000000e+01   2.82000008e+01 \
         3.97999992e+01   5.62000008e+01   7.94000015e+01 \
         1.12000000e+02   1.58000000e+02   2.24000000e+02 \
         3.16000000e+02   4.47000000e+02   6.31000000e+02 \
         8.91000000e+02   1.26000000e+03   1.78000000e+03 \
         2.51000000e+03   3.55000000e+03   5.01000000e+03 \
         7.08000000e+03   1.00000000e+04
do
	j=0
	for o in 25.79999924  36.90000153  45.59999847  53.09999847 60. 66.40000153 72.5  78.5 84.30000305
	do
		#[ -e $2_${i}_${j}_rdata.hdf5 ] ||
		python3 torusBN.py --nh=$nh --opening-angle=$o --nevents $1 --output="$2_${i}_${j}_" &
		((j++))
	done
	wait
	((i++))
done

exit

# if $2 == output/bntorus
cd output
python3 ../xspecexport/createtorustable.py bntorus.fits bntorus_*_?_rdata.hdf5
python3 ../xspecexport/createtorustable.py bntorus-transmit.fits bntorus_*_?_transmitrdata.hdf5
python3 ../xspecexport/createtorustable.py bntorus-reflect.fits bntorus_*_?_reflectrdata.hdf5

python3 ../xspecexport/createtoruscutofftable.py bntorus-cutoff.fits bntorus_*_?_rdata.hdf5 
python3 ../xspecexport/createtoruscutofftable.py bntorus-cutoff-transmit.fits bntorus_*_?_transmitrdata.hdf5
python3 ../xspecexport/createtoruscutofftable.py bntorus-cutoff-reflect.fits bntorus_*_?_reflectrdata.hdf5

