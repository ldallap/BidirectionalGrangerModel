# BidirectionalGrangerModel

This repository is part of the paper: Feedforward and feedback influences through distinct frequency bands between two spiking-neuron networks, by: Leonardo Dalla Porta, Daniel M. Castro, Mauro Copelli, Pedro V. Carelli, and Fernanda S Matias.

Physical Review E: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.054404

arXiv: https://arxiv.org/abs/2106.05783


Whenever this code is used, please cite the aforomentioned paper.

########

The code is implemeted in c++ language.
In order to generate the simulated LFP (Local Field Potential) and the 10-s long spike-trains run the code *BidirectionalGrangerModel.cpp*

The code can be compiled as:
*g++ BidirectionalGrangerModel.cpp -O3 -march=native -o exec*

In order to execute, run:
*./exec 123456*

*123456* is the random number seed, you can use whatever integer number here.

The code will generate two outputs (text file):
1) FeedBackFeedForward_MeanMemPot_PARAMETERDESCRIPTION.mp 
2) SpikeTrain_FeedBackFeedForward_PARAMETERDESCRIPTION.st

File #1 is composed of three columns: 1st) time (ms); 2nd) LFP from Pop2 (alpha oscillations); 3rd) LFP from Pop1 (gamma oscillations)
File #2 is composed of 200 columns where each column is a binary spike train for different neurons. Columns from 1-100 stands for neurons in Pop2 and from 101-200 for neurons in Pop1. Each line corresponds to one step of simulation.

Tested in: gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0

Contact: leonardodallaporta@gmail.com
