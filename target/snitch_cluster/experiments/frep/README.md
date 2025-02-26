Run RTL experiments:
```
./experiments.py none.yaml sw run visual-trace -j
```

Run PLS experiments:
```
make clean-vsim
make PL_SIM=1 DEBUG=ON VCD_DUMP=1 bin/snitch_cluster.vsim
./experiments.py pls.yaml run power -j --run-dir pls
```
