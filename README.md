# Parallel RePlAce

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

Disclaimer: reported performance is obtained from source, performance in container not measured.

### Docker

In `./`

`docker build . -t <your_tag>`

Then

`docker -it run <your_tag>`

Or 

`docker -it run <your_tag> <numberOfThreads> <benchmark_name> <csv_output>`

*Image size: ~10GB*

### Installed



## Experimental results

|          | RePlAce |         | Sequential |         | 1 Thread |         | 2 Threads |         | 4 Threads |         | 8 Threads |         | 12 Threads |         | RePlAce (12 t.) |         | DREAMPlace (12 t.) |         |
|----------|---------|---------|------------|---------|----------|---------|-----------|---------|-----------|---------|-----------|---------|------------|---------|-----------------|---------|--------------------|---------|
|          | GP      | TOTAL   | GP         | TOTAL   | GP       | TOTAL   | GP        | TOTAL   | GP        | TOTAL   | GP        | TOTAL   | GP         | TOTAL   | GP              | TOTAL   | GP                 | TOTAL   |
| adaptec1 | 108.68  | 189.30  | 59.60      | 140.21  | 58.68    | 137.97  | 47.57     | 117.91  | 39.43     | 103.64  | 35.59     | 94.64   | 35.22      | 93.33   | 46.45           | 104.07  | 171.13             | 198.63  |
| adaptec2 | 210.07  | 309.29  | 134.40     | 234.44  | 134.16   | 232.65  | 113.82    | 202.51  | 93.54     | 172.15  | 83.54     | 157.52  | 82.77      | 154.99  | 107.52          | 179.21  | 275.01             | 308.98  |
| adaptec3 | 443.67  | 628.37  | 294.82     | 481.08  | 298.33   | 484.75  | 231.02    | 395.15  | 179.49    | 328.32  | 151.36    | 290.37  | 144.50     | 282.42  | 215.93          | 352.01  | 457.26             | 520.85  |
| adaptec4 | 497.98  | 691.84  | 346.26     | 537.55  | 350.52   | 539.44  | 268.59    | 441.89  | 202.76    | 358.27  | 165.36    | 312.01  | 156.04     | 299.48  | 270.93          | 413.08  | 515.88             | 588.46  |
| bigblue1 | 205.09  | 308.17  | 106.40     | 209.47  | 107.72   | 209.15  | 84.17     | 174.03  | 69.56     | 150.98  | 62.83     | 137.58  | 60.88      | 133.85  | 89.25           | 161.91  | 238.53             | 275.14  |
| bigblue2 | 376.93  | 607.11  | 257.62     | 490.28  | 258.59   | 487.34  | 202.43    | 409.36  | 154.68    | 347.26  | 132.42    | 314.88  | 125.73     | 306.18  | 199.54          | 378.20  | 438.05             | 535.23  |
| bigblue4 | 2461.81 | 3478.51 | 1571.31    | 2580.64 | 1675.59  | 2678.92 | 1335.57   | 2212.92 | 1003.35   | 1781.11 | 821.90    | 1548.15 | 763.90     | 1472.50 | 1073.69         | 1777.88 | 2264.48            | 2610.96 |
| geomean  | 370.59  | 556.67  | 231.12     | 418.62  | 234.01   | 419.65  | 185.80    | 351.07  | 146.51    | 295.34  | 126.16    | 264.50  | 121.21     | 256.89  | 177.66          | 313.14  | 427.91             | 494.52  |


## Publication



## Citation

```bibtex

```

## Authors
- @gessfred Frédéric Gessler(frederic.gessler@epfl.ch) EPFL
- @mirjanastojilovic Mirjana Stojilovic(mirjana.stojilovic@epfl.ch) EPFL
- Philip Brisk UCSD