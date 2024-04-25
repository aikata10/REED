## CHIEF

This artifact is submitted to support our work: CHIEF: Chiplet-based Accelerator for Fully Homomorphic Encryption.

## Getting started

This is the root folder and contains the following files and directories. 


```
.
|__CKKS_CA_Model
|__NTT_INTT_unit
|__README.md
```

## CONTENT DESCRIPTION

1. Folder-`CKKS_CA_Model` contains our sample cycle accurate model used to estimate the performance. It features the possibility to run multiple chiplet configurations as multiple threads and a C2C (thread-to-thread) communication. We also show HBM data distribution along with task distribution for CKKS macro routines. The model is coupled with a working CKKS library to ensure correctness.
2. Folder- `NTT_INTT_unit` contains the NTT/INTT implementation, which closely follows Algorithm-1 in the paper. Along with the actual implementation, it also includes a toy example to help understand the data and control flow with proper visualization.



For more details on how to run these artefacts, please refer to the README.md file in the respective folders.


## License

CHIEF is licensed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3. See [LICENSE](LICENSE) for more details.



