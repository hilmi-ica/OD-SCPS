# Online Discovery and Decision Support for Subsea Seawater Injection Systems
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

## Reference

M. Hilmi, J. Vatn, A. Hasan, and M. A. Lundteigen, “Online model discovery under parameter uncertainty for risk and reliability-based decision support,” in Proceedings of the 36th European Safety and Reliability Conference, 2026.
```
@INPROCEEDINGS{OD-SCPS,
    author = {Hilmi, Muhammad and Vatn, Jørn and Hasan, Agus and Lundteigen, Mary Ann},
    title = {Online Model Discovery under Parameter Uncertainty for Risk and Reliability-Based Decision Support},
    booktitle = {Proceedings of the 36th European Safety and Reliability Conference},
    year = {2026}
}
```
