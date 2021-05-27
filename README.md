# BidirectionalGrangerModel

This repository is part of the paper: XXX

The code is implemeted in c++ language.
In order to generate the simulated LFP (Local Field Potential) and the 10-s long spike-trains run the code *BidirectionalGrangerModel.cpp*

The code can be compiled as:
*g++ BidirectionalGrangerModel.cpp -O3 -march=native -o exec*

To run the code:
*./exec Parameter1* (Parameter1 is a integer number which is used as the random number generator seed)

The code will generate two outputs (text file):
1) FeedBackFeedForward_MeanMemPot_PARAMETERDESCRIPTION.mp 
2) SpikeTrain_FeedBackFeedForward_PARAMETERDESCRIPTION.st

File #1 is composed of three columns: 1st) time (50000 ms); 2nd) LFP from Pop2 (alpha oscillations); 3rd) LFP from Pop1 (gamma oscillations)
File #2 is composed of 200 columns where each column is a binary spike train for different neurons. Columns from 1-100 stands for neurons in Pop2 and from 101-200 for neurons in Pop1. Each line corresponds to one step of simulation.


Contact: leonardodallaporta@gmail.com
