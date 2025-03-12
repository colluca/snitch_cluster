Run RTL experiments:
```
make clean-vsim
make bin/snitch_cluster.vsim DEBUG=ON -j
./experiments.py --actions sw run visual-trace -j
```

PLS simulations:

To test the PLS results
```
make clean-vsim
make PL_SIM=1 DEBUG=ON bin/snitch_cluster.vsim
./experiments.py --actions run -j --run-dir pls_test
```

To run the power simulation
```
make clean-vsim
make PL_SIM=1 DEBUG=ON VCD_DUMP=1 bin/snitch_cluster.vsim
./experiments.py power.yaml --actions run power -j --run-dir pls_power
```

To build all hardware configurations:
```
BIN_DIR=$PWD/experiments/frep/hw/base32fc/bin   make VSIM_BUILDDIR=experiments/frep/hw/base32fc/work-vsim   DEBUG=ON CFG_OVERRIDE=cfg/base32fc.hjson   $BIN_DIR/snitch_cluster.vsim
BIN_DIR=$PWD/experiments/frep/hw/zonl32fc/bin   make VSIM_BUILDDIR=experiments/frep/hw/zonl32fc/work-vsim   DEBUG=ON CFG_OVERRIDE=cfg/zonl32fc.hjson   $BIN_DIR/snitch_cluster.vsim
BIN_DIR=$PWD/experiments/frep/hw/zonl64fc/bin   make VSIM_BUILDDIR=experiments/frep/hw/zonl64fc/work-vsim   DEBUG=ON CFG_OVERRIDE=cfg/zonl64fc.hjson   $BIN_DIR/snitch_cluster.vsim
BIN_DIR=$PWD/experiments/frep/hw/zonl64dobu/bin make VSIM_BUILDDIR=experiments/frep/hw/zonl64dobu/work-vsim DEBUG=ON CFG_OVERRIDE=cfg/zonl64dobu.hjson $BIN_DIR/snitch_cluster.vsim
BIN_DIR=$PWD/experiments/frep/hw/zonl48dobu/bin make VSIM_BUILDDIR=experiments/frep/hw/zonl48dobu/work-vsim DEBUG=ON CFG_OVERRIDE=cfg/zonl48dobu.hjson $BIN_DIR/snitch_cluster.vsim
```

Command to run a single experiment (from the experiment's run directory) for debugging:
```
../../../../../../../../../sw/blas/gemm/scripts/verify.py ../../../../../hw/zonl48dobu/bin/snitch_cluster.vsim ../../../../../build/zonl48dobu/24/32/8/gemm.elf
```
