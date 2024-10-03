# Project Overview
This project was developed as part of a 3-year degree program at the Università degli Studi di Milano (University of Milan). It aims to provide a comprehensive implementation and analysis of the GGH lattice-based cryptosystem and its GGH-HNF variant using modern tools. The work explores the resilience of lattice-based cryptography against quantum threats, offering in-depth mathematical analysis, cryptanalysis, and experimental results on performance and security. A novel hybrid variant is introduced, combining elements of both systems. By comparing findings with classical cryptosystems and previous studies, the project assesses the practical viability of these lattice-based schemes in current technological contexts.
# GGH_crypto Python Package
This repository contains a Python package implementing the Goldreich-Goldwasser-Halevi (GGH) public key cryptosystem and its optimization, GGH-HNF by Micciancio. The GGH cryptosystem is a lattice-based cryptographic system proposed in 1997, which offers potential resistance against quantum computer attacks. Unfortunately, despite their theoretical interest, both the original GGH cryptosystem and its GGH-HNF optimization have been declared as practically unusable due to security vulnerabilities and performance limitations. This implementation serves primarily for educational and research purposes, offering insights into the challenges and developments in lattice-based cryptography.
The package uses and implements essential functions for lattice-based cryptography, such as:
- Algorithms for solving the Closest Vector Problem (CVP):
    - Babai's rounding technique
    - Embedding technique
- Lattice reduction algorithms:
    - Wrapped FPLLL implementation of BKZ
    - Gram-Schmidt orthogonalization
    - FLINT implementation of LLL

The package is composed by 3 specific subpackages:
- `GGH.py` which implements the original cryptosystem presented in 1997.
- `GGH_HNF.py` which implements the optimised version of GGH proposed by Micciancio in 2002.
- `Utils.py` which includes functions and methods in common between the 2 cryptosystems. It also contains general algorithms -or wrappers for them- related to the lattice-based cryptography, such as: BKZ, Gram-Schmidt, and CVP solvers.  

# Installation
The Package can be simply installed with PIP: 
```
pip install GGH-crypto
```
However, to use the BKZ reduction, it's necessary to install the external FPLLL library. FPLLL is only compatible with Linux. You can refer to [this link](https://github.com/fplll/fplll#compilation) for installation steps on various operative systems.

# General Usage
To demonstrate the usage and capabilities of this package, several example scripts are provided in the Examples folder:
1. `GGH_classes_example.py`: This script showcases the basic usage of the GGH and GGH-HNF classes, demonstrating key generation, encryption, and decryption processes.
2. `GGH_HNF_attack.py`: This example implements an attack on the GGH-HNF variant described in Christoph Ludwig's technical report presented in 2004.
3. `GGH_attack_Nguyen.py`: This script demonstrates the Nguyen attack on the original GGH cryptosystem, published in 1999. This example uses the simplest case where the public matrix is invertible mod 2σ.
4. `embedding_technique.py`: This example showcases the embedding technique for solving the Closest Vector Problem (CVP). It also includes a visualization of the problem using matplotlib.
5. `rounding_off_technique.py`: This script implements Babai's rounding-off algorithm, another approach to solving the CVP. Like the embedding technique example, it also provides a visualization of the problem using matplotlib.

Each function and its parameters are thoroughly described with comments in the package. Furthermore vectors within this package must adhere to the [FLINT Python bindings](https://fredrikj.net/python-flint/) notation for proper definition and representation.
The hybrid version proposed in the thesis, which outperforms both GGH and GGH-HNF in various fields, can be tested by setting the `GGH_private` flag to `True`. Doing so will generate a random private basis using the method applied in GGH, while following the process implemented in GGH-HNF for subsequent steps.

# References and final notes

Testing demonstrated significant performance improvements. The hybrid version showed superior results in both performance and security, requiring twice the time to break compared to the original GGH. However, with a dimension of 800 needed for security, these systems still present practical limitations compared to RSA, mainly in decryption times and key sizes. Nevertheless, the proposed hybrid version could open new perspectives for future analyses.

This thesis is sadly written in Italian only. For non-Italian speakers interested in reading the details, please use a translation tool to access the content. 

References can be found at the end of **thesis.pdf** file.
