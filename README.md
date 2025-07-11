# About
MWC Apple Apple Silicon is designed and optimized for Apple Silicon CPU, but in theory should work on amy device that support Metal.

# Build setup (Mac OS)

Install cmake if you don't have it:
```
> brew install cmake
```


Install jsoncpp:
```
> brew install jsoncpp
```


# Clone and build the miner
Note, mwc-metal-miner requires metal-cpp in the same dev root directory. That is why please clone both those projects
```
> cd <your_dev_root_directory>
> git clone https://github.com/bkaradzic/metal-cpp.git
> git clone https://github.com/mwcproject/mwc-metal-miner.git
> cd mwc-metal-miner
> cmake -DBUILD_TARGET=CPU_MINER -DCMAKE_BUILD_TYPE=Release .
> make
```

mwc-metal-miner executable should be located in your current directory. Run it to see the usage.

```
> ./mwc-metal-miner
Usage: ./mwc-metal-miner -node <host:port> -login <user_name> [-pass <password>]
```

# Run miner

If you are mining solo, please check documentaiton how to setup the node and the wallet: https://github.com/mwcproject/mwc-node/tree/master/doc/mining

If you are using mining pool, you need to know host:port for stratum connection and your user name.

cpu-mwc-miner usage:
```
Usage: ./mwc-metal-miner -node <host:port> -login <user_name> [-pass <password>]
```

-node:  host:port of mining node stratum interface. If you run your mining node locally with activated stratum protocol, use parameter `-node 127.0.0.1:3416` for mainnet and `-node 127.0.0.1:13416` for floonet.
-login: your user name for stratum or mining pool. If you are mining solo, you can use any name.
-pass: optional password for your login. Might be required by some mining pools.

