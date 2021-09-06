# MPA-OpenCl
### OpenCL MPA API ###  
===========================================================================

DESCRIPTION:

MPA-OpenCL is a Multiple precision Arithmetic API in OpenCl licensed under the Apache Software License v2.
The API offer a set of helper functions that can be used to carry usual Arbitrary Precision Arithmetic in every openCL enable device like GPUs, Multi core CPUs, Co-processeurs, FPGA and hand held devices that support OpenCL like for example Android Devices supporting OpenCL.
The main motivation behind the development of this API come from the lack of such an API in OpenCL, similar API exists for proprietary GPGPU platform like CUDA.
The API can be used to accelerate applications using multiple precision arithmetics like ECDSA, RSA, Research in physics, Big Data analysis Applications to name a few.

List of supported multiple precision arithmetic operations in the MPA-OpenCl version 1.0-beta

 1-Big numbers comparaison 
 
 2-Big numbers addition
 
 3-Big numbers subtraction 
 
 4-Big numbers multiplication (several algorithms)
 
 5-Big numbers exponentiation
 
 6-Big numbers division
 
 7-Big numbers integer square root
 
 7-Big numbers reduction
 
 8-Big numbers Modular addition
 
 9-Big numbers Modular subtraction 
 
10-Big numbers Modular multiplication 

11-Montgomery multiplication for big numbers

12-Big numbers Modular exponentiation

# Kernels Depensancies : 
MPA-OpenCL kernels use stansards OpenCL-C language as described in Khronos group specifications and do not use any external dependencies.

# Host Dependancies : 
The Host application also use standard C 99 and do not depends on external libraries 

# Test Dependencies : 
To test the application result and compare theme with other MPA API you will need the OpenSSL and GMP librairies.



