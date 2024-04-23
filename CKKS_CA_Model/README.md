
## CHIEF

Welcome to CHIEF Artefacts!

In support of our work: CHIEF: Chiplet-based Accelerator for Fully Homomorphic Encryption.

## Getting started

This directory contains the sample proof-of-work for the CKKS cycle-accurate model utilized in this work. It is coupled with a levelled CKKS implementation to ensure correct computation and runtime. 

Its contents are as follows.


```
.
|__Makefile
|__context_class.cpp 
|__main.cpp
|__Scheme.cpp
|__numb.h/cpp
|__data_types.cpp 
|__polynomial_functions_class.
|__number_theory_functions_class.cpp
|__/sage_generated_keys
|__README.md
```
 

## Prerequisite

The implementation requires a GNU C++ compiler for the cpp code and a Sage for generating keys (if new keys are needed).


## HOW TO RUN
To run, the user should execute the following command-

> $ make clean ; make ; ./main ; make clean


## CONTENTS

This folder contains several files and two main subfolders.

1. The subfolder- `sage_generated_keys` contains the pre-generated secret, public, and key-switching (evaluation) keys.
2. File- `context_class.cpp` contains the global parameters such as the ring degree, chiplet configuration, and depth L.
3. File- `main.cpp` is the base file that performs addition and multiplication with relinearization. It also reports runtime.
4. File- `Scheme.cpp` contains the multi-threaded instantiations of basic CKKS routines- Multiply, Add, KeySwitch, ModDown, and Rescale.
5. The remaining files contain polynomial-wise and coefficient-wise helper functions.


## Information on the parameters and chiplet-based design simulation

1. The chiplet behaviour is simulated using multi-threading routine openmp. 
2. The number of chiplets can be varied {1,2,4,8} by changing the variable `Num_Chiplets` in `context_class.cpp`.
3. Based on both the SW runtime and the cycle-accurate HW runtime, the performance benefit of chiplet-based FHE HW design is validated.

## Information on HBM data distribution and task distribution

The evaluation keys are distributed across chiplets at the very beginning of the computation. This is shown in file `data_types.cpp` in function 'HBM(Context& context, Key& key, unsigned seed)'. 

This comprises the initial data distribution. After this, ciphertexts to be processed are received by the chiplets distributed limb-wise (as detailed in the paper). This is then processed parallelly by all the chiplets for all the different operations offered in the CKKS scheme (`Scheme.cpp`).

We have shown the parallel chiplet-based task distribution for macro operations- Add, Mult, KeySwitch, ModDown, and ReScale in `Scheme.cpp`. Each operation is accompanied by cycle-accurate modelling and computation verification to ensure the cycles are obtained for correct execution only.

