// for the `memset`
#include <string.h>
// generating sine and cosine samples
#include <math.h>

#include "fft.h"



// build and run:		(with Linux - also link with math: -lm)
//		gcc example.c -o example && ./example



static void print_input_time_domain_samples(float complex in[], size_t n, const char *title) {
	printf("\n\t%s\n", title);
	for (size_t i = 0UL; i < n; i++) {
		printf("%zu ", i);

		// presumed known amplitude - so that negative parts can be visible and printed
		float t = in[i] * 10.f + 30.f;
		for (int i = 0; i < t - 1; i++) printf(" ");
		printf("*");
		for (int i = 45 - 1; i > t; i--) printf(" ");
		printf("\n");
	}
}

static void print_frequency_domain_samples(float complex fft_out[], size_t n) {
	printf("\n\tFFT MAGNITUDE\n");
	const float fft_factor = 10.f;		// used so that stems are visible after normalization (to scale them up)

	for (size_t i = 0UL; i < n; i++) {
		printf("%zu ", i);
		float a = fft_factor * 2.f * cabsf(fft_out[i]) / n;
		if (a >= 0.001f) { for (int i = 0; i < a; i++) printf("x"); }
		printf("\n");
	}
	printf("\n");



	printf("\n\tFFT PHASE\n");				// in radians
	for (size_t i = 0UL; i < n; i++) {
		printf("# ");
		float a = cabsf(fft_out[i]);
		float phi = fft_factor * cargf(fft_out[i]);
		if (a >= 0.001f) { for (int i = 0; i < a; i++) printf("x"); }
		printf("\n");
	}
	printf("\n");



	// complex form of fft samples
	// for (int i = 0; i < n; i++) printf("%2.2f + %2.2fi  ", crealf(fft_out[i]), cimagf(fft_out[i])); printf("\n");
}



int main() {
	size_t n = 1<<6;
	float complex in[n];
	float complex fft_out[n];
	float complex ifft_out[n];

	memset(fft_out, 0, n);
	memset(ifft_out, 0, n);

	for (size_t i = 0; i < n; i++) {
		float t = (float) i / n;
		in[i] = 2.f * cosf(2*PI*t*1) + sinf(2*PI*t*4);
	}



	// original samples
	print_input_time_domain_samples(in, n, "ORIGINAL SAMPLES");

	// first invoke requires stride to be 1 (further internal calls use it)
	fft(in, 1, fft_out, n);

	// frequency spectrum
	// print_frequency_domain_samples(fft_out, n);

	
	// inverse fft samples
	ifft(fft_out, 1, ifft_out, n);
	for (size_t i = 0UL; i < n; i++) ifft_out[i] *= 1.f / (float) n;		// normalization

	print_input_time_domain_samples(ifft_out, n, "RECOVERED SAMPLES");


	// // converting complex to real - sort of signed cabsf
	// float ifft_real[n];
	// for (size_t i = 0UL; i < n; i++) ifft_real[i] = ifft_out[i];
	

	return 0;
}