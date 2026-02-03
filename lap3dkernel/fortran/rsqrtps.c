#include <immintrin.h>

// Approximate reciprocal sqrt for 4 floats with two Newton refinements in double.
// out[i] ~= 1/sqrt(in[i])
void rsqrtps_nr4(const float in[4], double out[4]) {
    __m128 x = _mm_loadu_ps(in);
    __m128 y = _mm_rsqrt_ps(x);
    float y_f[4];
    _mm_storeu_ps(y_f, y);
    for (int i = 0; i < 4; ++i) {
        double xd = (double)in[i];
        double yd = (double)y_f[i];
        // Newton-Raphson iterations in double:
        // y = y * (1.5 - 0.5*x*y*y)
        yd = yd * (1.5 - 0.5 * xd * yd * yd);
        yd = yd * (1.5 - 0.5 * xd * yd * yd);
        out[i] = yd;
    }
}

// Approximate reciprocal sqrt for 1024 floats using rsqrtps + two Newton refinements in double.
void rsqrtps_block1024(const float *in, double *out) {
    const __m256d half = _mm256_set1_pd(0.5);
    const __m256d three_halves = _mm256_set1_pd(1.5);
    for (int i = 0; i < 1024; i += 4) {
        __m128 x = _mm_loadu_ps(in + i);
        __m128 y = _mm_rsqrt_ps(x);
        __m256d x_d = _mm256_cvtps_pd(x);
        __m256d y_d = _mm256_cvtps_pd(y);
        __m256d y2 = _mm256_mul_pd(y_d, y_d);
        y_d = _mm256_mul_pd(y_d, _mm256_sub_pd(three_halves, _mm256_mul_pd(_mm256_mul_pd(x_d, y2), half)));
        y2 = _mm256_mul_pd(y_d, y_d);
        y_d = _mm256_mul_pd(y_d, _mm256_sub_pd(three_halves, _mm256_mul_pd(_mm256_mul_pd(x_d, y2), half)));
        _mm256_storeu_pd(out + i, y_d);
    }
}
