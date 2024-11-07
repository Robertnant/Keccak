# SHA-3 Hash Implementation

This project implements the SHA-3 cryptographic hash function based on the Keccak algorithm. It allows for hashing files with various SHA-3 hash sizes, including 224, 256, 384, and 512 bits.

## Features

- **Keccak Permutation**: Implements the theta, rho, pi, chi, and iota steps for the Keccak-f permutation.
- **Sponge Construction**: Utilizes the sponge function for absorbing input and squeezing out the desired hash output.
- **Configurable Hash Size**: Supports SHA-3 hash sizes of 224, 256, 384, and 512 bits.
- **File Hashing**: Can compute the SHA-3 hash of specified files.

## Requirements

- Python 3
- NumPy library (`pip install numpy`)

## Usage

To run the script, use the following command:

```bash
python SHA3.py -m [224|256|384|512] file1 [file2 ...]
```

### Arguments

- `-m [224|256|384|512]`: Specify the SHA-3 hash size.
- `file1 [file2 ...]`: Paths to the files to be hashed.

### Example

```bash
python SHA3.py -m 256 example.txt
```

This command will compute the SHA3-256 hash of `example.txt`.

## Code Overview

- **theta, rho, pi, chi, iota**: These functions implement the core steps of the Keccak-f permutation.
- **keccak_p and keccak_f**: Functions to apply Keccak-f with the specified number of rounds.
- **sponge**: Sponge construction for hashing, using `keccak_f` for compression.
- **pad10star1**: Padding function to align data to the block size.
- **sha3**: Computes the SHA-3 hash of a message.
- **main**: Parses command-line arguments and computes the hash for specified files.
