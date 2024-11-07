from collections.abc import Sequence
from math import ceil
from math import log2
import sys
import os
import numpy as np

# 1st step of the Keccak-f permutation
def theta(State):
  w = 1600 // 25

  # Compute the parity 
  C = [
      [
          State[x][0][z] ^ State[x][1][z] ^ State[x][2][z] ^ State[x][3][z] ^ State[x][4][z]
          for z in range(len(State[0][0]))
      ]
      for x in range(len(State))
  ]

  # XOR the parities of the previous and the rotated next lane
  D = [
      [
          C[(x - 1) % 5][z] ^ C[(x + 1) % 5][(z - 1) % w]
          for z in range(len(State[0][0]))
      ]
      for x in range(len(State))
  ]

  # XOR the original state with the XORed parities
  State_prime = [
      [
          [State[x][y][z] ^ D[x][z] for z in range(len(State[0][0]))]
          for y in range(len(State[0]))
      ]
      for x in range(len(State))
  ]
  return State_prime

#Usual values for the rotation offsets for each round
rotation_values = [1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300]

def rho(State):
    w = 1600 // 25

    # Initialize State_prime with all zeroes
    State_prime = [[[0 for _ in range(w)] for _ in range(5)] for _ in range(5)]

    # Copy the first lane
    State_prime[0][0] = State[0][0][:]
    
    # Coordinates (x, y) for each step
    coordinates = [(1, 0)]
    for _ in range(1, 24):
        x, y = coordinates[-1]
        coordinates.append((y, (2 * x + 3 * y) % 5))

    # Apply rotations
    for i, (x, y) in enumerate(coordinates):
        rotation = rotation_values[i]
        State_prime[x][y] = [State[x][y][(z - rotation) % w] for z in range(w)]

    return State_prime

# Coordinates for the pi step
pi_coordonees = [(0, 0), (3, 0), (1, 0), (4, 0), (2, 0), (1, 1), (4, 1), (2, 1), (0, 1), (3, 1), (2, 2), (0, 2), (3, 2), (1, 2), (4, 2), (3, 3), (1, 3), (4, 3), (2, 3), (0, 3), (4, 4), (2, 4), (0, 4), (3, 4), (1, 4)]

def pi(State):
    size = len(State[0][0])
    State_prime = [[[0 for _ in range(size)] for _ in range(5)] for _ in range(5)]

    # Apply the permutation
    for i in range(25):
            new_x, new_y = pi_coordonees[i]
            for z in range(size):
                State_prime[i//5][i%5][z] = State[new_x][new_y][z]
    return State_prime

def chi(State):
  # Apply the chi step to the state array
  State_prime = [
      [
          [
              State[x][y][z] ^ ((State[(x + 1) % 5][y][z] ^ 1) & State[(x + 2) % 5][y][z])
              for z in range(len(State[0][0]))
          ]
          for y in range(len(State[0]))
      ]
      for x in range(len(State))
  ]
  return State_prime

# Round constants for the iota step
round_constants = [
    0x0000000000000001, 0x0000000000008082, 0x800000000000808A,
    0x8000000080008000, 0x000000000000808B, 0x0000000080000001,
    0x8000000080008081, 0x8000000000008009, 0x000000000000008A,
    0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
    0x000000008000808B, 0x800000000000008B, 0x8000000000008089,
    0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
    0x000000000000800A, 0x800000008000000A, 0x8000000080008081,
    0x8000000000008080, 0x0000000080000001, 0x8000000080008008
]

def iota(State, ir):
    w = 1600 // 25
    round_constant = round_constants[ir]
    # Applying the round constant to the first lane of the state
    for z in range(w):
        State[0][0][z] ^= (round_constant >> z) & 1
    return State

def keccak_p(b, num_rounds, S):
  w = b // 25

  # 1. Let State = (S[0], S[1], …, S[24]), where S[i] is a 64-bit string.
  State = [
      [[S[w * (5 * y + x) + z] for z in range(w)] for y in range(5)]
      for x in range(5)
  ]
  # 2. For ir from 12+2l –nr to 12 + 2l –1, let State = Rnd(State, ir).
  l = int(ceil(log2(w)))  # Default should be 6.

  # 3. For ir from 0 to 23, let State = ι(χ(π(ρ(θ(State)))), ir).  
  for ir in range(12 + 2 * l - num_rounds, 12 + 2 * l):
    State = iota(chi(pi(rho(theta(State)))), ir)

  # 4. Assign value to new state
  A_out = np.zeros(b, dtype=int)  
  for x in range(5):
    for y in range(5):
      for z in range(w):
        A_out[w * (5 * y + x) + z] = State[x][y][z]
  return A_out 

def keccak_f(b, S):
  # Launch 24 rounds of the Keccak-f permutation
  return keccak_p(b, 24, S)

def sponge(f, pad, r, N, d):
  # 1. Pad the input
  P = np.concatenate((N, pad(r, len(N))))
  n = len(P) // r
  c = 1600 - r  # b - r
  P_seq = np.split(P, n)
  w = 1600 // 25  # b / 25
  S = np.zeros(1600, dtype=int)

  # 2. Absorb phase
  for i in range(n):
    Pi = P_seq[i]
    s = np.bitwise_xor(S, np.concatenate([Pi, np.zeros(c, dtype=int)]))
    S = f(1600, s)

  # 3. Squeeze phase
  Z = S[:r]
  while d > len(Z):
    S = f(1600, S)
    Z = np.concatenate([Z, S[:r]])
  return Z[:d]

def pad10star1(x, m):
  # 1. Compute pad string
  j = (-m - 2) % x
  return np.array([1] + [0] * j + [1])

def keccak(c, N, d):
  # Launch the sponge function
  return sponge(keccak_f, pad10star1, 1600 - c, N, d)

def reverse_bits_in_bytes(bit_array):
  # Reverse the bits in each byte of the array
  reversed_bytes = bit_array.reshape(-1, 8)[:, ::-1].flatten()
  return reversed_bytes

def sha3(message, hashsize):
  # Compute the SHA3 hash of the message
  binary_result = keccak(hashsize*2, np.concatenate((message, np.array([0, 1]))), hashsize)
  binary_result = reverse_bits_in_bytes(binary_result)
  reversed_hex = np.packbits(binary_result).tobytes().hex()
  return reversed_hex

def main(argv):
  # Check if the first argument is '-m'
  if len(argv) < 2 or argv[1] != '-m':
    print("Error: First argument must be '-m'")
    print("Usage: python script.py -m [224|256|384|512] file1 [file2 ...]")
    return
  # Get sha3 version to use
  hash_size = 256  # Default
  if len(argv) > 2:
    if argv[2] in ["224", "256", "384", "512"]:
        hash_size = int(argv[2])
    else:
        print("Invalid SHA3 size. Using default 256.")

  # Compute the hash of each file
  for file_path in argv[3:]:
    if not os.path.exists(file_path):
      print(f"File not found: {file_path}")
      continue

    with open(file_path, "rb") as file:
      file_content = open(file_path, "rb").read()
      if len(file_content) == 0:
          bit_array = np.array([], dtype=int)
      else:
          bit_array = np.array([(b >> i) & 1 for b in file_content for i in range(8)])
      result = sha3(bit_array, hash_size)
      print(f"SHA3-{hash_size} ({file_path}) = {result}")

if __name__ == "__main__":
  main(sys.argv)
