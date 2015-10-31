FacSexBGSSims:

Forward-in-time simulations of background selection, that can account for varying rates 
of selfing and asexuality. Written in C++.

This simulation is based on those created for Kamran-Disfani & Agrawal 2014 JEB paper
"Selfing, adaptation and background selection in finite populations".

Variables are changed in the file 'theConstants.h'; the program is then compiled and executed. 
By default, a burn-in population is created for 2N generations, then the population is sampled at 200 timepoints over a
further 2N generations. Various outputs are produced, including variance in the neutral
quantitative trait (which is used to measure Ne).

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.
