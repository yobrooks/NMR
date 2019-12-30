#Nuclear Magnetic Resonance Spectrum

##Table of Contents
* [Introduction](#introduction)
* [Features](#features)
* [Setup and Execution](#setup-and-execution)

##Introduction
This project analyzes discrete data from Nuclear Magnetic Resonance (NMR) spectrometer. The purpose of this project is to fit the given data to a natural cubic spline (after putting the discrete points through a filtering technique to smooth the data and remove noise from it) and then find the peak locations in the spline. This program was developed using C++ for the final project component of a Numerical Methods class in Fall 2019.

##Features
To run the program needs an input file which will hold all of the user specified options for the NMR analysis. The default name of the input file is "nmr.in". To change how the program analyzes the NMR spectrum, change the input options in the "nmr.in" input file. The program also needs a ".dat" file to hold the discrete data points that make up the spectrum. The default is "testdata2.dat".

The program has four different options for filtering the discrete data: no filter, boxcar, Savitsky-Golay, and Discrete Fourier. The data can be passed through the specified filter any number of times as the user wishes. The Discrete Fourier Transform filtering technique uses Linear Algebra and has three different solvers to recover the filtered y points: Inverse Solver, Direct Solver (LU Decomposition), and Iterative Solver (Jacobi Iterative Method).

The program has four different options for integration techniques to find the areas under of the peaks: Adaptive Quadrature, Romberg, Simpson, and Gaussian Quadrature.

The program also finds the Tetramethylsilane (TMS) peak which is the most positive x-axis peak in the spectrum. Once this value is found, all other values are adjusted relative to it so that the maximum of the TMS peak can occur at 0.0 on the x-axis of the spectrum.


##Setup and Execution
To compile this program in a Unix/Linux environment, type `make` in the console. Then type `./main` in the console to run it. 

Upon the end of the program, the user is prompted to view whatever file was entered in as the
output file name to see the analysis and results of the NMR process.

In the default case it is "analysis.txt". To view the output file type `cat output_file_name`.
