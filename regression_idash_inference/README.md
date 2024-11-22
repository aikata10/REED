## Logistic Regression over SEAL

# Build
Run the following
```console
foo@bar:~$ cmake -S . -B build
foo@bar:~$ cmake -- build build
```

# Run
```console
foo@bar:~$ cd build
foo@bar:~$ make 
foo@bar:~$ ./regression_clean 
```

# main.cpp

This is the only source code of this project, and it requires seal.

# Dataset

The dataset is given in idash.txt and is from idash 2017.

The first column is cancer status, where 1 means positive and 0 means not.

Other columns are used as inputs of regression, you can check the results by copy - of course, except for the first data (cancer status) - and paste it into vData in line 60 of main.cpp.
