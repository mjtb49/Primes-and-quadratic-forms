# Primes of the form x^2 + ny^2
A criteria to identify precisely the primes of the form x^2 + ny^2 is described in a book by Cox. This repository computes that criteria.

The main difficult thing this program must do is compute the minimal polynomial for a primitive element of the Hilbert class field lying over an imaginary quadratic extension. For this purpose, the theory of complex multiplication allows us to build abelian extensions of such fields, and the Hilbert class field, as the maximal unramified abelian extension, turns out the be generated the values of the Klein J invariant on representatives of the ideal classes. Therefore, we are able to approximate the minimal polynomial by computing the Klein J invariant on each representative, then finding a polynomial with integer coefficients having these as roots. The fact that this relies on numerical computation of the J function is unfortunate, but in practice the computed values for this repository tend to be remarkably close to integers for n < 2000, and one can update the precision found in the main method if one needs bigger n

Worth noting is that for some numbers, the 65 "convenient numbers" of which the largest is 1848, the polynomial is not necessary and a mere condition on the remainder of the prime modulo the discriminant suffices.

In principle a small modification to this code should be able to determine conditions for primes of the form x^2 + xy + ny^2
