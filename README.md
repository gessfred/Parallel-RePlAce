# Parallel RePlAce

Based on commit `dbbe12c` from [The-OpenROAD-Project](https://github.com/The-OpenROAD-Project/RePlAce).

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Based on the paper "A Shared-Memory Parallel Implementation of the RePlAce Global Cell Placer" from VLSID2020.

RePlAce is a state-of-the-art prototype of a flat, analytic, and nonlinear global cell placement algorithm, which models a placement instance as an electrostatic system with positively charged objects. 
This extension includes techniques to reduce memory contention and to effectively balance the workload among threads, targeting the most substantial performance bottlenecks. 
With 2–12 threads, our parallel RePlAce speeds up the bin density function by a factor of 4.2–10×, the wirelength function by a factor of 2.3–3×, and the cost gradient function by a factor of 2.9–6.6× compared to the single-threaded original RePlAce baseline. Moreover, our parallel RePlAce is ≈3.5× faster than the state-of-the-art PyTorch-based placer DREAMPlace, when both are running on 12 CPU cores

## Installation

Install the dependencies:

- flex
- bison
- Boost
- Intel IPP
- Intel MKL
- X11

For Debian:

```shell
apt-get install -y libstdc++6 wget build-essential libx11-dev libboost-dev libcurl4 cmake swig flex bison libtool zlib1g-dev tcl-dev tk-dev libjpeg-dev
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
sh -c 'echo deb https://apt.repos.intel.com/ipp all main > /etc/apt/sources.list.d/intel-ipp.list'
apt-get update -y
apt-get install -y intel-mkl-2018.2-046 intel-ipp-2018.4-057 
```

Then go to `code/RePlAce` and do `make`

## Usage

### Docker

In `./`

`docker build . -t <your_tag>`

Then

`docker -it run <your_tag>`

Or 

`docker -it run <your_tag> <numberOfThreads> <benchmark_name> <csv_output>`

*Image size: ~10GB*

### Installed



## Publication

[preprint](https://infoscience.epfl.ch/record/273833?ln=en)

## Citation

```bibtex

```

## Authors
- @gessfred Frédéric Gessler(frederic.gessler@epfl.ch) EPFL
- @mirjanastojilovic Mirjana Stojilovic(mirjana.stojilovic@epfl.ch) EPFL
- Philip Brisk UCSD