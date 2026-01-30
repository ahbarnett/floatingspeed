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
