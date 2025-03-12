Run RTL experiments:
```
make clean-vsim
make bin/snitch_cluster.vsim DEBUG=ON -j
./experiments.py --actions hw sw run visual-trace -j
```

PLS simulations:

To test the PLS results:
```
make clean-vsim
make PL_SIM=1 DEBUG=ON bin/snitch_cluster.vsim
./experiments.py --actions run -j --run-dir pls_test
```

To run the power simulation:
```
make clean-vsim
make PL_SIM=1 DEBUG=ON VCD_DUMP=1 bin/snitch_cluster.vsim
./experiments.py power.yaml --actions run power -j --run-dir pls_power
```

Command to run a single experiment (from the experiment's run directory) for debugging:
```
../../../../../../../../../sw/blas/gemm/scripts/verify.py ../../../../../hw/zonl48dobu/bin/snitch_cluster.vsim ../../../../../build/zonl48dobu/24/32/8/gemm.elf
```
