# Online Discovery and Decision Supprt for Subsea Systems
This repository contains simulations for submersible cascade pump system for online identification and decision support tools using sparse WyNDA (with recursive sparse regression and recursive sparse Bayesian learning), supporting our submission to the European Safety and Reliability Conference.

## Purpose

Sharing MATLAB code and data for reproducing results and facilitating review.

## Repository Structure

* **`ODWyNDA.m`**: MATLAB simulation code for online discovery of the sparse WyNDA framework.
  * **`WILibFuncs.m`**: Basis library functions to accompany `ODWyNDA.m`.
  * **`NWIData.mat`**: Sequential simulation data to accompany `ODWyNDA.m`.
* **`DSWyNDA.m`**: MATLAB simulation code of the sparse WyNDA framework for decision support applications.
  * **`BBNFuncs.m`**: Bayesian belief network function to accompany `DSWyNDA.m`.
  * **`FWIData.mat`**: Sequential simulation data to accompany `DSWyNDA.m`.

## European Safety and Reliability Conference

This repository directly supports data and code for our paper submission.
