# DLCA 3D

## Overview
This code written as part of my master's project.

It simulates the stochastic process of aggregation between particles.

## About the code files

Use make to compile the C++ source files. It will generate 2 executables - simulation & multirun.
Usage:

```bash
./simulation <size> <con1> <con2> <runs> <sticking>
```

size - Lattice size (Simulation space size)
con1 - concentration of particle type 1
con2 - concentration of particle type 2
runs - number of simulation runs (default = 1)
sticking - probability that particles of type 1 and 2 stick.



Use the python files to plot graphs and visualize.

## Example

```bash
make
./simulation 10 0.064

Container Size = 10
particle concentration = 0.064
polymer concentration = 0
Number of runs = 1
Sticking Probability = 1


Simulation 1
Number of molecules is 64

Starting Simulation with 53 clusters

Simulation 1 Done
All simulations done
```

Move the plotnano.py file to *System_States/system1* and run it. You'll get a visualization like folowwing.

![64 particles - 10x10](https://user-images.githubusercontent.com/666808/52552879-93d40200-2d96-11e9-8d70-770908887978.png)

There are a bunch of other files generated in system folders.
nano.txt - co-ordinates of particles of type 1
plmr.txt - co-ordinates of particles of type 2

*I'll try to add more documentation later if possible*