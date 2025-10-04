# GeoHAR-SPDNet
Bucci, Palma, Zhang - GeoHAR SPDNet Version 1.0,  Copyright(c) November, 2024. 

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [License](#license)


## Installation
To get started, clone this repository and install the required dependencies.

```bash
git clone https://github.com/C4ntor/GeoHAR-SPDNet.git
cd GeoHAR-SPDNet
```

Ensure Machine Learning toolbox is installed

## Usage
Step1: Place your data in vech form under the folder "./data". Note that the data used in the paper "RCOVReal.csv" are already provided in data folder.

Step2: Launch spdnet_afew.m for a simple example. Results willbe saved in the root folder of the project, or launch neuralnet_afew.m to get benchmarks for geometric-less predictions for the network that is not compliant to the manifold geometry.  For this latter it is also possible to remove regeig by editing neuralnet_init_afew.m as follows:  

```bash
opts.layernum = 1; 
...
net.layers = {} ;
net.layers{end+1} = struct('type', 'bfc',...
                          'weight', Winit{1}) ;
%net.layers{end+1} = struct('type', 'rec') ;
%net.layers{end+1} = struct('type', 'bfc',...
%                          'weight', Winit{2}) ;
%net.layers{end+1} = struct('type', 'rec') ;
%net.layers{end+1} = struct('type', 'bfc',...
%                          'weight', Winit{3}) ;
net.layers{end+1} = struct('type', opts.loss_function) ;



## License
Note that the copyrights of the manopt toolbox, and SPDNet are reserved respectively by Manopt <a href="https://www.manopt.org/"></a> and Zhiwu Huang and Luc Van Gool. A Riemannian Network for SPD Matrix Learning, In Proc. AAAI 2017 <a href="https://github.com/zhiwu-huang/SPDNet"></a>

## How to Cite <a name="How-to-Cite"></a>





