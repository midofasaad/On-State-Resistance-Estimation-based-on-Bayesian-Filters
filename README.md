# On-State-Resistance-Estimation-based-on-Bayesian-Filters
The core goal of this project implements the kalman filter and  particle filter for the estimation of the sinusodial current and voltage signal.  

### Table of Contents

1. [Instructions](#instructions)
2. [Project Motivation](#motivation)
3. [Licensing, Authors, and Acknowledgements](#licensing)

## Instructions: <a name="instructions"></a>

The project was created with Python 3.9.0.
The following files has the following functions:

1.  openCmd.bat creates and activate virtual environment in folder env38.
2. 'Kalman Filters.py' implments and visualize an estimation from the kalman filter for a sinusodial signal. The intial and structural parameters of the filter and the signal can be modified.
by the user of the code. 
3. 'Particle Filters.py' implments and visualize an estimation from the particle filter for  a sinusodial signal. The intial and structural parameters of the filter and the signal can be modified.
by the user of the code. 
4. 'Parameter Estimation.py' visualize a normal Distribution with arbitrary mean and variance, such that parameter estimation methods (e.g maximum liklihood) can be tested to determine
             the parameters of the distribution.

## Project Motivation: <a name="motivation"></a>

An accurate estimation of the current and voltage measurement signals allows the accurate estimation of the on-state resistance. This allows the usage of the developed lookup tables
for the real time estimation of junction temperature. 


## Licensing, Authors, Acknowledgements: <a name="licensing"></a>

This corporate message data is from one of the free datasets provided on the
[Figure Eight Platform](https://appen.com/resources/datasets/), licensed under
a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

Feel free to use my code as you please:

Copyright 2020 Mahmoud Saad

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
