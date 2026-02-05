// madpcm.h
// A more modern approach to ADPCM.
#pragma once
#ifndef MADPCM_H
#define MADPCM_H

#include <stdbool.h>
#include <stdint.h>

#define MADPCM_LPC_ORDER 4
#define MADPCM_BLOCK_SAMPLES 1024
#define MADPCM_HEADER_SIZE 24
#define MADPCM_PAYLOAD_SIZE 512
#define MADPCM_MONO_SIZE (MADPCM_HEADER_SIZE + MADPCM_PAYLOAD_SIZE) 
#define MADPCM_STEREO_SIZE (MADPCM_MONO_SIZE * 2)
#define MADPCM_WAVE_FORMAT 0x4D41 // 'MA' in hex

// Core MADPCM functions
int madpcm_encode_block(const int16_t* inSamples, const int16_t* prevSamples, uint8_t* outBuffer, bool fast);
int madpcm_encode_block_stereo(const int16_t* inL, const int16_t* inR, const int16_t* prevL, const int16_t* prevR, uint8_t* outBuffer, bool fast);
void madpcm_decode_block(const uint8_t* inBuffer, int16_t* outSamples);
void madpcm_decode_block_stereo(const uint8_t* inBuffer, int16_t* outL, int16_t* outR);

#ifdef MADPCM_IMPLEMENTATION

// prepare for the clusterfuck
#pragma region freestanding
#ifndef MADPCM_FREESTANDING
#include <string.h>
#define madpcm_memset memset
#define madpcm_memcpy memcpy
#else // freestanding
// i am running on a toaster implementation
#ifndef MADPCM_MEMFUNCS
void* memcpy(void* dest, const void* src, size_t count);
void* memset(void* dest, int val, size_t count);
#define madpcm_memset memset
#define madpcm_memcpy memcpy
#else
// truly freestanding
#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(__movsb, __stosb)
#endif
static inline void* madpcm_memcpy(void* dest, const void* src, size_t count) {
#if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
    __movsb((unsigned char*)dest, (const unsigned char*)src, count);
    return dest;
#elif defined(__GNUC__)
    return __builtin_memcpy(dest, src, count);
#else
    uint8_t* d = (uint8_t*)dest; const uint8_t* s = (const uint8_t*)src;
    while (count--) *d++ = *s++;
    return dest;
#endif
}
static inline void* madpcm_memset(void* dest, int val, size_t count) {
#if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
    __stosb((unsigned char*)dest, (unsigned char)val, count);
    return dest;
#elif defined(__GNUC__)
    return __builtin_memset(dest, val, count);
#else
    uint8_t* d = (uint8_t*)dest; while (count--) *d++ = (uint8_t)val;
    return dest;
#endif
}
#endif
#endif
#pragma endregion

#define MADPCM_TRELLIS_K 4
#define MADPCM_Q15_MAX 32767
#define MADPCM_Q15_MIN -32768

// Gaussian Table (Q11) - Optimized for Laplacian residuals
static const int16_t madpcm__gauss_steps[8] = { 205, 614, 1024, 1433, 1843, 2252, 2662, 3200 };

typedef struct madpcm__block_header {
    int16_t scales[4]; // LSB of scales[0] is Stereo Flag
    int16_t coeffs[MADPCM_LPC_ORDER];
    int16_t seed[MADPCM_LPC_ORDER]; // History samples
} madpcm__block_header;

// Trellis helper struct
typedef struct madpcm__trellis_candidate {
    int64_t err;
    int32_t p_idx;
    int16_t val;
    int8_t  nib;
} madpcm__trellis_candidate;

typedef struct madpcm__trellis_state {
    int64_t err;
    int16_t h[MADPCM_LPC_ORDER];
} madpcm__trellis_state;

// Math Helpers
static inline int16_t madpcm__clamp_q15(int32_t x) {
    if (x > MADPCM_Q15_MAX) return MADPCM_Q15_MAX;
    if (x < MADPCM_Q15_MIN) return MADPCM_Q15_MIN;
    return (int16_t)x;
}
static inline int32_t madpcm__abs_int32(int32_t x) { return (x < 0) ? -x : x; }
static inline int64_t madpcm__abs_int64(int64_t x) { return (x < 0) ? -x : x; }

static void madpcm__compute_autocorr(const int16_t* signal, int len, int order, int64_t* r) {
    for (int i = 0; i <= order; i++) {
        int64_t sum = 0;
        for (int j = 0; j < len - i; j++) sum += (int64_t)signal[j] * signal[j + i];
        r[i] = sum;
    }
    r[0] = (r[0] > 256) ? r[0] : 256;
    r[0] += (r[0] >> 9);
}

static int madpcm__compute_LPC(const int64_t* r, int order, int16_t* outCoeffs) {
    int32_t a[MADPCM_LPC_ORDER + 1][MADPCM_LPC_ORDER + 1] = { 0 };
    int64_t e[MADPCM_LPC_ORDER + 1];
    e[0] = r[0];

    for (int k = 1; k <= order; k++) {
        int64_t sum = 0;
        for (int j = 1; j < k; j++) {
            int64_t prod = (int64_t)a[k - 1][j] * r[k - j];
            sum += (prod >> 24);
        }
        int64_t num = r[k] - sum;
        int64_t den = e[k - 1];
        if (den == 0) den = 1;
        if (madpcm__abs_int64(num) >= madpcm__abs_int64(den)) return 0;

        int32_t lambda = (int32_t)((num << 24) / den);
        a[k][k] = lambda;

        for (int j = 1; j < k; j++) {
            int64_t term = (int64_t)lambda * a[k - 1][k - j];
            a[k][j] = a[k - 1][j] - (int32_t)(term >> 24);
        }

        int64_t l2 = (int64_t)lambda * lambda;
        e[k] = (e[k - 1] * ((1LL << 30) - (l2 >> 18))) >> 30;
    }

    for (int i = 0; i < order; i++) {
        outCoeffs[i] = madpcm__clamp_q15(a[order][i + 1] >> 12);
    }
    return 1;
}

// LUT-based Dequantizer
static inline int16_t madpcm__dequantize_LUT(int8_t nibble, int16_t scale) {
    int8_t idx = nibble & 0x07;
    int8_t sign = (nibble & 0x08) ? -1 : 1;
    int32_t val = ((int32_t)scale * madpcm__gauss_steps[idx]) >> 11;
    return (sign < 0) ? (int16_t)-val : (int16_t)val;
}

// Finds best nibble in O(1) using inverse quantization estimation
static inline int8_t madpcm__find_best_nibble(int32_t target, int16_t scale) {
    if (scale < 1) scale = 1;

    int32_t absTarget = madpcm__abs_int32(target);
    // Approximation: index ~ (val * 160) >> 16
    // 1/410 (step size) in Q16 is approx 160
    int32_t val = (absTarget << 11) / scale;
    int32_t idx = (val * 160) >> 16;

    if (idx > 7) idx = 7;

    // Check neighbor to handle boundary cases
    int32_t est = (scale * madpcm__gauss_steps[idx]) >> 11;
    int32_t diff1 = madpcm__abs_int32(absTarget - est);

    if (idx < 7) {
        int32_t est2 = (scale * madpcm__gauss_steps[idx + 1]) >> 11;
        int32_t diff2 = madpcm__abs_int32(absTarget - est2);
        if (diff2 < diff1) idx++;
    }

    if (target < 0) return (int8_t)(idx | 0x08);
    return (int8_t)idx;
}

static inline void madpcm__stereo_to_mid_side(int16_t* l, int16_t* r, int count) {
    for (int i = 0; i < count; i++) {
        int32_t L = l[i], R = r[i];
        l[i] = (int16_t)((L + R) >> 1);
        r[i] = (int16_t)((L - R) >> 1);
    }
}

static inline void madpcm__mid_side_to_stereo(int16_t* l, int16_t* r, int count) {
    for (int i = 0; i < count; i++) {
        int32_t M = l[i], S = r[i];
        l[i] = madpcm__clamp_q15(M + S);
        r[i] = madpcm__clamp_q15(M - S);
    }
}

// Encoder
int madpcm_encode_block(const int16_t* inSamples, const int16_t* prevSamples, uint8_t* outBuffer, bool fast) {
    madpcm__block_header* header = (madpcm__block_header*)outBuffer;
    uint8_t* payload = outBuffer + sizeof(madpcm__block_header);

    // Capture History from Raw PCM
    // internal history layout: [0]=t-1, [1]=t-2, [2]=t-3, [3]=t-4
    if (prevSamples) {
        header->seed[0] = prevSamples[3]; 
        header->seed[1] = prevSamples[2];
        header->seed[2] = prevSamples[1]; 
        header->seed[3] = prevSamples[0];
    } else {
        madpcm_memset(header->seed, 0, sizeof(header->seed));
    }

    // High Precision Analysis
    int64_t r[MADPCM_LPC_ORDER + 1] = { 0 };
    madpcm__compute_autocorr(inSamples, MADPCM_BLOCK_SAMPLES, MADPCM_LPC_ORDER, r);

    if (!madpcm__compute_LPC(r, MADPCM_LPC_ORDER, header->coeffs)) {
        madpcm_memset(header->coeffs, 0, sizeof(header->coeffs));
    }

    // Simulation & Heuristics Loop
    int64_t sigEnergy = 0, resEnergy = 0;
    int zeroCrossings = 0;
    int32_t maxResSub[4] = { 0 };
    int32_t maxSigSub[4] = { 0 };
    
    // Initialize simulation history with perfect seeds
    int16_t simHist[MADPCM_LPC_ORDER];
    madpcm_memcpy(simHist, header->seed, sizeof(simHist));

    int32_t c0 = header->coeffs[0], c1 = header->coeffs[1];
    int32_t c2 = header->coeffs[2], c3 = header->coeffs[3];
    int badModel = 0;

    for (int i = 0; i < MADPCM_BLOCK_SAMPLES; i++) {
        int subIdx = i >> 8;

        // ZCR Check: (current ^ prev) < 0 means sign bit changed
        if (i > 0 && ((inSamples[i] ^ inSamples[i - 1]) < 0))
            zeroCrossings++;

        int32_t pred = (c0 * simHist[0]) + (c1 * simHist[1]) + (c2 * simHist[2]) + (c3 * simHist[3]);
        pred >>= 12;

        int32_t res = inSamples[i] - pred;
        sigEnergy += (int64_t)inSamples[i] * inSamples[i];
        resEnergy += (int64_t)res * res;

        if (madpcm__abs_int32(res) > maxResSub[subIdx]) maxResSub[subIdx] = madpcm__abs_int32(res);
        if (madpcm__abs_int32(inSamples[i]) > maxSigSub[subIdx]) maxSigSub[subIdx] = madpcm__abs_int32(inSamples[i]);

        if (madpcm__abs_int32(res) > 28000) badModel = 1;

        simHist[3] = simHist[2]; simHist[2] = simHist[1]; simHist[1] = simHist[0]; simHist[0] = inSamples[i];
    }

    int highZCR = zeroCrossings > (MADPCM_BLOCK_SAMPLES >> 2);
    int lowGain = sigEnergy < (8 * (resEnergy + 1));
    int useFallback = badModel || (highZCR && lowGain);

    if (useFallback) {
        madpcm_memset(header->coeffs, 0, sizeof(header->coeffs));
        c0 = c1 = c2 = c3 = 0;
    }

    // Populate scales with Gaussian Correction (~0.66x)
    for (int k = 0; k < 4; k++) {
        int32_t targetMax = useFallback ? maxSigSub[k] : maxResSub[k];
        int32_t s = (targetMax * 2) / 3;
        if (s < 10) s = 10;
        header->scales[k] = madpcm__clamp_q15(s);
    }

    // Trellis Quantization
    int16_t curHist[MADPCM_LPC_ORDER];
    madpcm_memcpy(curHist, header->seed, sizeof(curHist));

    // Process 4 sub-blocks
    for (int sb = 0; sb < 4; sb++) {
        int16_t scale = header->scales[sb];
        const int16_t* blk = inSamples + (sb * 256);
        uint8_t* subPayload = payload + (sb * 128);

        // FAST MODE: No Trellis
        if (fast) {
            for (int i = 0; i < 256; i += 2) {
                // Nibble 1
                int32_t p1 = (c0 * curHist[0] + c1 * curHist[1] + c2 * curHist[2] + c3 * curHist[3]) >> 12;
                int8_t n1 = madpcm__find_best_nibble(blk[i] - p1, scale);
                int16_t r1 = madpcm__clamp_q15(p1 + madpcm__dequantize_LUT(n1, scale));
                curHist[3] = curHist[2]; curHist[2] = curHist[1]; curHist[1] = curHist[0]; curHist[0] = r1;

                // Nibble 2
                int32_t p2 = (c0 * curHist[0] + c1 * curHist[1] + c2 * curHist[2] + c3 * curHist[3]) >> 12;
                int8_t n2 = madpcm__find_best_nibble(blk[i + 1] - p2, scale);
                int16_t r2 = madpcm__clamp_q15(p2 + madpcm__dequantize_LUT(n2, scale));
                curHist[3] = curHist[2]; curHist[2] = curHist[1]; curHist[1] = curHist[0]; curHist[0] = r2;

                subPayload[i >> 1] = (n1 & 0x0F) | (n2 << 4);
            }
            continue;
        }

        // QUALITY MODE: Trellis
        // Traceback buffers [sample][path]
        int8_t trace_nib[256][MADPCM_TRELLIS_K];
        int8_t trace_idx[256][MADPCM_TRELLIS_K];

        madpcm__trellis_state states[MADPCM_TRELLIS_K];
        int activePaths = 1;
        states[0].err = 0;
        madpcm_memcpy(states[0].h, curHist, sizeof(curHist));

        for (int i = 0; i < 256; i++) {
            madpcm__trellis_candidate cands[MADPCM_TRELLIS_K * 3];
            int candCount = 0;

            for (int p = 0; p < activePaths; p++) {
                int16_t* h = states[p].h;
                int32_t pred = (c0 * h[0] + c1 * h[1] + c2 * h[2] + c3 * h[3]) >> 12;
                int32_t target = blk[i] - pred;

                // Check best greedy + neighbors
                int8_t bestN = madpcm__find_best_nibble(target, scale);
                int8_t centerIdx = (bestN & 0x08) ? (bestN & 0x07) + 8 : bestN;

                for (int offset = -1; offset <= 1; offset++) {
                    int idx = centerIdx + offset;
                    if (idx < 0 || idx > 15) continue;

                    int8_t nib = (idx < 8) ? idx : (idx & 0x07) | 0x08;
                    int16_t rec = madpcm__clamp_q15(pred + madpcm__dequantize_LUT(nib, scale));
                    int64_t diff = blk[i] - rec;
                    int64_t sqErr = diff * diff;

                    cands[candCount].err = states[p].err + sqErr;
                    cands[candCount].p_idx = p;
                    cands[candCount].nib = nib;
                    cands[candCount].val = rec;
                    candCount++;
                }
            }

            // Select best K
            madpcm__trellis_state nextStates[MADPCM_TRELLIS_K];
            int nextCount = 0;

            for (int k = 0; k < MADPCM_TRELLIS_K && k < candCount; k++) {
                int bestJ = -1;
                int64_t minE = -1;
                for (int j = 0; j < candCount; j++) {
                    if (cands[j].err == -1) continue;
                    if (bestJ == -1 || cands[j].err < minE) { minE = cands[j].err; bestJ = j; }
                }
                if (bestJ == -1) break;

                // Save traceback data instead of full path copy
                trace_idx[i][k] = (int8_t)cands[bestJ].p_idx;
                trace_nib[i][k] = cands[bestJ].nib;

                nextStates[k] = states[cands[bestJ].p_idx];
                nextStates[k].err = cands[bestJ].err;

                int16_t* h = nextStates[k].h;
                h[3] = h[2]; h[2] = h[1]; h[1] = h[0]; h[0] = cands[bestJ].val;

                cands[bestJ].err = -1; // Mark used
                nextCount++;
            }
            activePaths = nextCount;
            madpcm_memcpy(states, nextStates, sizeof(madpcm__trellis_state) * activePaths);
        }

        // Commit best path (Backtracking)
        madpcm_memcpy(curHist, states[0].h, sizeof(curHist));
        int bestPathIdx = 0; // Index 0 is best due to sort

        uint8_t nibbles[256];
        for (int i = 255; i >= 0; i--) {
            nibbles[i] = trace_nib[i][bestPathIdx];
            bestPathIdx = trace_idx[i][bestPathIdx];
        }

        for (int j = 0; j < 128; j++) {
            subPayload[j] = (nibbles[2 * j] & 0x0F) | (nibbles[2 * j + 1] << 4);
        }
    }

    return sizeof(madpcm__block_header) + (MADPCM_BLOCK_SAMPLES / 2);
}

// Adaptive Stereo Encoder
int madpcm_encode_block_stereo(const int16_t* inL, const int16_t* inR, const int16_t* prevL, const int16_t* prevR, uint8_t* outBuffer, bool fast) {
    int16_t bufL[MADPCM_BLOCK_SAMPLES], bufR[MADPCM_BLOCK_SAMPLES];
    int16_t histL[4], histR[4];

    madpcm_memcpy(bufL, inL, sizeof(bufL));
    madpcm_memcpy(bufR, inR, sizeof(bufR));

    // Prepare local history copies (safe to modify)
    if (prevL && prevR) {
        madpcm_memcpy(histL, prevL, 4 * sizeof(int16_t));
        madpcm_memcpy(histR, prevR, 4 * sizeof(int16_t));
    } else {
        madpcm_memset(histL, 0, sizeof(histL));
        madpcm_memset(histR, 0, sizeof(histR));
    }

    // Adaptive Check: Calculate Energies
    int64_t eLR = 0, eMS = 0;
    for (int i = 0; i < MADPCM_BLOCK_SAMPLES; i += 4) {
        int32_t m = (bufL[i] + bufR[i]) >> 1;
        int32_t s = (bufL[i] - bufR[i]) >> 1;
        eLR += madpcm__abs_int32(bufL[i]) + madpcm__abs_int32(bufR[i]);
        eMS += madpcm__abs_int32(m) + madpcm__abs_int32(s);
    }

    bool useMS = (eMS < eLR * 0.9f);
    
    // If using MS, we must also convert the HISTORY to MS
    if (useMS) {
        madpcm__stereo_to_mid_side(bufL, bufR, MADPCM_BLOCK_SAMPLES);
        madpcm__stereo_to_mid_side(histL, histR, 4); 
    }

    int s1 = madpcm_encode_block(bufL, histL, outBuffer, fast);
    int s2 = madpcm_encode_block(bufR, histR, outBuffer + s1, fast);

    // Flagging: Use LSB of first scale
    madpcm__block_header* h = (madpcm__block_header*)outBuffer;
    if (useMS) h->scales[0] |= 1;
    else       h->scales[0] &= ~1;

    return s1 + s2;
}

// Decoder
void madpcm_decode_block(const uint8_t* inBuffer, int16_t* outSamples) {
    const madpcm__block_header* header = (const madpcm__block_header*)inBuffer;
    const uint8_t* payload = inBuffer + sizeof(madpcm__block_header);

    const int32_t c0 = header->coeffs[0];
    const int32_t c1 = header->coeffs[1];
    const int32_t c2 = header->coeffs[2];
    const int32_t c3 = header->coeffs[3];

    int16_t h0 = header->seed[0];
    int16_t h1 = header->seed[1];
    int16_t h2 = header->seed[2];
    int16_t h3 = header->seed[3];

    int sample = 0;

    // Each scale applies to 256 samples
    for (int s = 0; s < MADPCM_BLOCK_SAMPLES / 256; s++) {
        int16_t scale = header->scales[s] & ~1;

        for (int i = 0; i < 256; i++, sample++) {
            int32_t pred =
                (c0 * h0 + c1 * h1 + c2 * h2 + c3 * h3) >> 12;

            uint8_t byte = payload[sample >> 1];
            int8_t nibble = (sample & 1)
                ? (byte >> 4)
                : (byte & 0x0F);

            int16_t delta = madpcm__dequantize_LUT(nibble, scale);
            int32_t val = pred + delta;

            h3 = h2;
            h2 = h1;
            h1 = h0;
            h0 = madpcm__clamp_q15(val);

            outSamples[sample] = h0;
        }
    }
}

void madpcm_decode_block_stereo(const uint8_t* inBuffer, int16_t* outL, int16_t* outR) {
    const madpcm__block_header* h = (const madpcm__block_header*)inBuffer;
    bool isMS = (h->scales[0] & 1);

    int sz = sizeof(madpcm__block_header) + (MADPCM_BLOCK_SAMPLES / 2);
    madpcm_decode_block(inBuffer, outL);
    madpcm_decode_block(inBuffer + sz, outR);

    if (isMS) madpcm__mid_side_to_stereo(outL, outR, MADPCM_BLOCK_SAMPLES);
}

#endif // MADPCM_IMPLEMENTATION
#endif // MADPCM_H