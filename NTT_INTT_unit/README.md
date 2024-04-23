## CHIEF

Welcome to CHIEF Artefacts!

In support of our work: CHIEF: Chiplet-based Accelerator for Fully Homomorphic Encryption.

## Getting started

This directory structure contains the files for the new design of the proposed novel *HYBRID* NTT/INTT unit. 


```
.
|__helper.py
|__hybrid_ntt_demo.py
|__parameters.py
|__README.md
```


# CONTENT DESCRIPTION
> File-`parameters.py` contains two sets of parameters. First is a toy parameter with a small polynomial degree (64) and prime modulus (257) to help understand the workings of our NTT/INTT unit. It can be disabled by setting the variable `USE_TOY_EXAMPLE` to false within this file. Doing this will set the polynomial degree to 2^16 (65536) and a bigger prime modulus (9007199256051713).

> File- `hybrid_ntt_demo.py` contains the NTT/INTT implementation, which very closely follows Algorithm 1 in the paper. For the toy example, proper prints are given to help visualize the data flow at every step. 

> File-`helper.py` contains the necessary smaller modules for NTT/INTT computation.


# HOW TO RUN

To run the NTT/INTT
> $ python3 hybrid_ntt_demo.py 

To change parameters, dive into parameters.py. The current parameters can be seen by running the following:
> $ python3 parameters.py 
