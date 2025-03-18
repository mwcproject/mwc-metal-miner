# About
MWC CPU miner is designed and optimized for Apple Silicon CPU. See details below about performance and memory consumtpion.

# Build setup (Mac OS)
Install jsoncpp:
```
brew install jsoncpp
```

# Compile test and version

## Tests 

Tests are expected to be run in the Debug mode. Only developers are expected to run them if code was updated.
First build the tests and run them.
```
cmake -DBUILD_TARGET=TESTS -DCMAKE_BUILD_TYPE=Debug .
make
./mwc_cpu_miner
```
Running tests will take about 15 minutes. Internal states are validated with 'assert' calls. 

## Build miner

```
cmake -DBUILD_TARGET=CPU_MINER -DCMAKE_BUILD_TYPE=Release .
make
./mwc_cpu_miner
```

# Run miner

cpu-mwc-miner usage:
```
./mwc_cpu_miner -node <host:port> -login <user_name> [-pass <password>] -algo <C31|C32>
```

'-node' should point to the mwc-node that runs with enable stratum. 
