Run RTL experiments:
```
make clean-vsim
make bin/snitch_cluster.vsim
./experiments.py --actions sw run visual-trace -j
```

PLS simulations:

To test the PLS results
```
make clean-vsim
make PL_SIM=1 DEBUG=ON bin/snitch_cluster.vsim
./experiments.py --actions run -j --run-dir pls_check
```

To run the power simulation
```
make clean-vsim
make PL_SIM=1 DEBUG=ON VCD_DUMP=1 bin/snitch_cluster.vsim
./experiments.py power.yaml --actions run power -j --run-dir pls
```
