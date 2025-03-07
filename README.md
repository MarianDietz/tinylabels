# TinyLabels

This artifact provides the implementation used to evaluate the efficiency claims of the paper [TinyLabels: How to Compress Garbled Circuit Input Labels, Efficiently](https://eprint.iacr.org/2024/2048), published at EUROCRYPT 2025.

**Do not use this code in a productive environment; no security guarantees are made.**

## Overview

The code is based on a patch of the Microsoft SEAL framework, which is used for operations on ring elements as required by RingLWE.
The original SEAL codebase has been modified to allow for plaintext modulus that divides is ring modulus.

This implementation does not support all settings described in the paper. It only supports moduli that are the product of two different NTT-friendly primes, where one of them is simultaneously as the plaintext modulus.
No RO-trick for compression of ct2 has been implemented.

See the evaluation section of the paper for more details on the precise parameter setting.

## Requirements & Setup

Building the source code follows the same procedure as the SEAL framework (see README-SEAL.md).
In particular, it is sufficient to run the following two commands to build the implementation in-place:
```
cmake -S . -B build -DSEAL_BUILD_TINYLABELS=ON
cmake --build build
```

This will generate all executables (`setup`, `enc1`, `enc2`, `keygen`, `dec`, `gen_samples`, `benchmark`) of our implementation in the `build/bin` directory.
To test our implementation, run `cd build/bin`, and then follow the next steps.

If your processor supports [Intel HEXL](https://github.com/intel/hexl), you may enable the acceleration library using `-DSEAL_USE_INTEL_HEXL=ON` when running the first `cmake` commands. Doing so may yield some performance improvements.

## How to Use

After building the code, all algorithms of the batch-select primitive can be executed individually.

We recommend to first run `./gen_samples`: This will generate uniformly random vectors `l_1` and `l_2`, as well as a uniformly random binary vector `y`.
They are placed in files `l_1.txt`, `l_2.txt`, and `y.txt`.
Furthermore, the "expected" outcome of running the entire batch-select pipeline (i.e., `l_1 * y + l_2`) is placed into the file `expected.txt`.
You may modify the input vectors, or choose them entirely by yourself instead of running `./gen_samples`.
However, in order for the remaining algorithms to run without errors, it is necessary that the input files contain exactly `2^21` numbers (unless hardcoded parameters have been changed).

Then, the following algorithms should be executed (in this order):
* `./setup`: This generates the public parameters, saved into `pp.bin`.
* `./enc1`: This encrypts the vector given in `l_1.txt`, and the output is saved into `ct1.bin` (note that `pp.bin` must have been generated already). Furthermore, a state `st1.bin` is created, needed for key generation later.
* `./enc2`: This encrypts the vector given in `l_2.txt`, and the output is saved into `ct2.bin` (note that `pp.bin` must have been generated already). Furthermore, a state `st2.bin` is created, needed for key generation later.
* `./keygen`: If files `pp.bin`, `st1.bin` and `st2.bin` exist, this executable takes the binary vector given in the file `y.txt`, and saves the key into `sk.bin`.
* `./dec`: If files `pp.bin`, `ct1.bin`, `ct2.bin`, and `sk.bin` exist, this executable decrypts the evaluated ciphertext, and saves it into `output.txt`.

In order to verify that execution was correct, you can compare `expected.txt` with the actual decryption output `output.txt`, for example by running `diff expected.txt output.txt`.

The algorithms above can be repeatedly executed, for example to use `enc2` for encrypting several different vectors `l_2`, which may be useful to amortize the cost of `ct1.txt` across several instances of batch-select as described in the paper.

All executables will output statistics. To verify the efficiency claims made in the paper, compare the output "Total time" with row "Time" of Table 4 in the ePrint paper.
Each algorithm also outputs its running time split into the different types of ring element operations. These values correspond to those listed in Table 5 in the ePrint paper.

Note that storage space is not optimized in our implementation (e.g., 128 bits are required to store a single polynomial coefficient of bitlength 109), and therefore the sizes of `ct1.bin` etc. are slightly larger than the sizes claimed in the paper.
The main purpose of this implementation is the evaluation of running time.

## Modifying Parameters

The constants at the beginning of `native/mytests/batchselect.h` may be modified to test the implementation on other parameters.
For example, by changing `w` to another value, the input vector length will be changed to `w*poly_modulus_degree`.
After re-building the source code and re-running `./gen_samples`, you can observe how the running time changes for the desired vector length.
