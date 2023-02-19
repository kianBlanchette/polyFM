# polyFM
This is a set of R functions used to estimate the carrier frequency of a sinusoidal signal modulated in frequency by a low-degree polynomial. They are slightly updated versions of the functions written during my time as a masters student at Queen's University. 

We take a semi-parametric approach to the signal detection, based on Thomson's multitaper method. The algorithm is, loosely, as follows:

-Given a time series with stationary noise. 

-Compute the eigencoefficients for the series.

-Project back to the time domain to create the "standard inverse" and its derivative.

-Use these to compute the instantaneous frequency series as an estimate of the modulating function.

-Compute the eigencoefficients of the instantaneous frequency series.

-Compute polynomial coefficients at degrees of your choice.

-Compute F statistics from the instantaneous frequency eigencoefficents and polynomial coefficients.

-Test the hypothesis that there is not a sinusoidal signal with frequency f modulated by a polynomial of degree p for all chosen degrees p and all frequencies f in the band of interest.

There are some added utility functions for computing the Harmonic F statistic and finding statistically significant local maxima of the test statistics, as well as some very untested and unfinished functions for jackknifing. I may return to explore the use of jackknifing some more eventually.

More context for what these functions can actually be used for can be found here: https://qspace.library.queensu.ca/handle/1974/28118 
or here: https://ieeexplore.ieee.org/abstract/document/9414209

Dependencies: the multitaper R package
