# FFT-related functions

One-file collection of few main FFT-related functions:
- `fft` - Fast Fourier Transform algorithm (translated from [Python](https://rosettacode.org/wiki/Fast_Fourier_transform#Python) to C)
- `ifft` - Inverse FFT (done with a neat [trick](https://stackoverflow.com/questions/17349424/is-it-possible-to-derive-a-the-2d-inverse-fft-algorithm-using-an-existing-1d-fft))
- `fft2` - 2D Fourier Transform (used with i.e. images)
- `ifft2` - Inverse FFT2
- `fft_shift` - shifts all rows and columns so that, in frequency domain, low frequencies move to the center of the 2D matrix (related to 2D transform - `fft2`)
- `fft_ishift` - undos `fft_shift`

Feel free to let me know if some fragments are unnecessary, optimizable or buggy. Had developed this for generating pink noise in mind.

## Using in projects

Straight up `#include "fft.h"`. Not sure if `<complex.h>` needs linking with anything (specifically with `-lm` on Linux).

## Notes

- FFT algorithm expects the `in` and `out` arrays to be equal sizes. Moreover, since it works best with sizes that are powers of 2, the safest for outcomes is `memset`-ing `out` to 0, so that you don't have to encounter garbage from memory, and deal with values like -2.43511e17.
- Results of `fft` and `ifft` are not normalized.
- Discrete time transform is periodic, with a period of sample rate. In `example.c` correct results appear for frequencies up to half the size of transformed array, then they are affected by results from the next iteration.