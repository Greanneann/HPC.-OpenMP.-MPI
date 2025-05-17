# Lattice Boltzmann Method (LBM) — MPI Communication Project

## Overview
This project implements communication routines for a parallel Lattice Boltzmann solver using MPI.

## Key Implementations

- `lbm_comm_init_ex1()`: 1D domain split; communicator size checks; rank assignment.
- `lbm_comm_ghost_exchange_ex1()`: Blocking left/right column exchange.
- `lbm_comm_ghost_exchange_ex2()`: Odd-even communication (even ranks receive first).
- `lbm_comm_ghost_exchange_ex3()`: Non-blocking communication with `MPI_Isend/Irecv` and `MPI_Waitall`.

- `lbm_comm_init_ex4()`: 
  - Introduces 2D domain splitting.
  - 2D Cartesian grid created with `MPI_Cart_create`.
  - Buffers allocated for up/down exchanges.
- `lbm_comm_release_ex4()`: Frees allocated buffers and communicator.
- `lbm_comm_ghost_exchange_ex4()`:
  - Same left/right exchange.
  - Manual buffer packing for up/down exchanges.
- `lbm_comm_init_ex5()`:
  - Adds derived MPI datatype (`MPI_Type_vector`) for vertical communication.
- `lbm_comm_ghost_exchange_ex5()`:
  - Uses derived types for up/down communication.
- `lbm_comm_ghost_exchange_ex6()`:
  - Non-blocking communication for 2D case (left/right and up/down exchanges).

## How to Run
- use `make`.
- Run with `mpirun -np 4 ./lbm -e 'number_of_exercise`.

