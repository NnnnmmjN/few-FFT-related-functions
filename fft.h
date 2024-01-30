#include <stdio.h>
#include <stdlib.h>
#include <complex.h>



#ifndef PI
#define PI 3.14159265359f
#endif



/*
// as the slower version - placed for reference
void dft(float complex in[], float complex out[], size_t n) {
	for (size_t f = 0; f < n; f++) {
		out[f] = 0;
		for (size_t i = 0; i < n; i++) {
			float t = (float) i / n;
			out[f] += in[i] * cexp(2*PI*f * t * I);
		}
	}
}
*/

void fft(float complex in[], size_t stride, float complex out[], size_t n) {
	if (n <= 1) { out[0] = in[0]; return; }

	fft(in,			 stride*2, out,		  n/2);
	fft(in + stride, stride*2, out + n/2, n/2);

	for (size_t k = 0; k < n/2; k++) {
		float t = (float) k / n;
		float complex v = out[k + n/2] * cexpf(-2.f*PI*I * t);
		float complex e = out[k];
		out[k] = e + v;
		out[k + n/2] = e - v;
	}
}

void ifft(float complex in[], size_t stride, float complex out[], size_t n) {
	(void) stride;
	float complex *temp = calloc(2*n, sizeof(float complex));
	float complex *temp_out = &temp[n];
	if (temp == NULL) return;

	for (size_t i = 0UL; i < n; i++) temp[i] = conjf(in[i]);
	fft(temp, 1, temp_out, n);
	for (size_t i = 0UL; i < n; i++) out[i] = conjf(temp_out[i]);

	if (temp) free(temp);
}

void fft2(float complex in[], float complex out[], size_t n, size_t data_stride) {
	float complex *rows_out = calloc(2*n, sizeof(float complex));
	float complex *cols_out = &rows_out[n];
	if (rows_out == NULL) return;

	size_t rows = n / data_stride;

	// row-wise fft
	for (size_t i = 0UL; i < rows; i++)
		fft(in + i*data_stride, 1, rows_out + i*data_stride, data_stride);
	// transpose
	for (size_t i = 0UL; i < n; i++)
		cols_out[i] = rows_out[(i % rows) * data_stride + (i / rows)];
	// column-wise fft (row-wise after transposition)
	for (size_t i = 0UL; i < data_stride; i++)
		fft(cols_out + i*rows, 1, rows_out + i*rows, rows);
	// transpose back to original
	for (size_t i = 0UL; i < n; i++)
		out[i] = rows_out[(i % data_stride) * rows + (i / data_stride)];

	if (rows_out) free(rows_out);
}

void ifft2(float complex in[], float complex out[], size_t n, size_t data_stride) {
	float complex *rows_out = calloc(2*n, sizeof(float complex));
	float complex *cols_out = &rows_out[n];
	if (rows_out == NULL) return;

	size_t rows = n / data_stride;
	// last transpose
	for (size_t i = 0UL; i < n; i++)
		rows_out[i] = in[(i % rows) * data_stride + (i / rows)];
	// column-wise ifft (row-wise after transposition)
	for (size_t i = 0UL; i < data_stride; i++)
		ifft(rows_out + i*rows, 1, cols_out + i*rows, rows);
	// transpose back
	for (size_t i = 0UL; i < n; i++)
		rows_out[i] = cols_out[(i % data_stride) * rows + (i / data_stride)];
	// row-wise ifft
	for (size_t i = 0UL; i < rows; i++)
		ifft(rows_out + i*data_stride, 1, out + i*data_stride, data_stride);

	if (rows_out) free(rows_out);
}

void fft_shift(float complex freqs[], size_t n, size_t stride) {
	float complex *tmp = calloc(n, sizeof(float complex));
	if (tmp == NULL) return;
	
	size_t rows = n / stride;
	size_t shift_x = stride / 2;
	size_t shift_y = rows / 2;

	for (size_t i = 0UL; i < n; i++) tmp[(i/stride) * stride + (shift_x + i) % stride] = freqs[i];
	for (size_t i = 0UL; i < n; i++) freqs[(shift_y * stride + i) % n] = tmp[i];

	if (tmp) free(tmp);
}

void fft_ishift(float complex freqs[], size_t n, size_t stride) {
	float complex *tmp = calloc(n, sizeof(float complex));
	if (tmp == NULL) return;
	
	size_t rows = n / stride;
	size_t shift_x = stride / 2;
	size_t shift_y = rows / 2;

	for (size_t i = 0UL; i < n; i++) tmp[i] = freqs[(shift_y * stride + i) % n];
	for (size_t i = 0UL; i < n; i++) freqs[i] = tmp[(i/stride) * stride + (shift_x + i) % stride];

	if (tmp) free(tmp);
}