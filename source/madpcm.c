// madpcm.c
// A more modern approach to ADPCM.
#include <stdbool.h>
#include <stdint.h>

#pragma region freestanding
#ifndef MADPCM_FREESTANDING
#include <string.h>
#include <stdlib.h>
#else // freestanding
// i am running on a toaster implementation
#ifndef MADPCM_MEMFUNCS
void* memcpy(void* dest, const void* src, size_t count);
void* memset(void* dest, int val, size_t count);
#else
// truly freestanding
#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(__movsb, __stosb)
#endif
static void* madpcm_memcpy(void* dest, const void* src, size_t count) {
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
static void* madpcm_memset(void* dest, int val, size_t count) {
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
#define memcpy madpcm_memcpy
#define memset madpcm_memset
#endif 
#endif 
#pragma endregion

#define LPC_ORDER 4
#define BLOCK_SAMPLES 1024
#define TRELLIS_K 4
#define Q15_MAX 32767
#define Q15_MIN -32768

// Gaussian Table (Q11) - Optimized for Laplacian residuals
static const int16_t kGaussSteps[8] = { 205, 614, 1024, 1433, 1843, 2252, 2662, 3200 };

typedef struct {
    int16_t scales[4]; // LSB of scales[0] is Stereo Flag
    int16_t coeffs[LPC_ORDER];
} MADPCMBlockHeader;

// Trellis helper struct
typedef struct {
    int64_t err;
    int pIdx;
    int8_t nib;
    int16_t val;
} TrellisCandidate;

typedef struct {
    int64_t err;
    int16_t h[LPC_ORDER];
    uint8_t path[256]; // Sub-block path
} TrellisPath;

// Math Helpers
static inline int16_t clamp_q15(int32_t x) {
    if (x > Q15_MAX) return Q15_MAX;
    if (x < Q15_MIN) return Q15_MIN;
    return (int16_t)x;
}

static inline int32_t abs_int32(int32_t x) { return (x < 0) ? -x : x; }
static inline int64_t abs_int64(int64_t x) { return (x < 0) ? -x : x; }

static void computeAutocorr(const int16_t* signal, int len, int order, int64_t* r) {
    for (int i = 0; i <= order; i++) {
        int64_t sum = 0;
        for (int j = 0; j < len - i; j++) sum += (int64_t)signal[j] * signal[j + i];
        r[i] = sum;
    }
    r[0] = (r[0] > 256) ? r[0] : 256;
    r[0] += (r[0] >> 9);
}

static int computeLPC(const int64_t* r, int order, int16_t* outCoeffs) {
    int32_t a[LPC_ORDER + 1][LPC_ORDER + 1] = { 0 };
    int64_t e[LPC_ORDER + 1];
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
        if (abs_int64(num) >= abs_int64(den)) return 0;

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
        outCoeffs[i] = clamp_q15(a[order][i + 1] >> 12);
    }
    return 1;
}

// LUT-based Dequantizer
static inline int16_t dequantizeLUT(int8_t nibble, int16_t scale) {
    int8_t idx = nibble & 0x07;
    int8_t sign = (nibble & 0x08) ? -1 : 1;
    int32_t val = ((int32_t)scale * kGaussSteps[idx]) >> 11;
    return (sign < 0) ? (int16_t)-val : (int16_t)val;
}

// Helper to find nearest Gaussian nibble
static int8_t findBestNibble(int32_t target, int16_t scale) {
    int32_t bestErr = 0x7FFFFFFF;
    int8_t bestN = 0;
    for (int i = 0; i < 16; i++) {
        int8_t packed = (i < 8) ? i : (i & 0x07) | 0x08;
        int16_t guess = dequantizeLUT(packed, scale);
        int32_t err = abs_int32(target - guess);
        if (err < bestErr) { bestErr = err; bestN = packed; }
    }
    return bestN;
}

static inline void stereoToMidSide(int16_t* l, int16_t* r, int count) {
    for (int i = 0; i < count; i++) {
        int32_t L = l[i], R = r[i];
        l[i] = (int16_t)((L + R) >> 1);
        r[i] = (int16_t)((L - R) >> 1);
    }
}

static inline void midSideToStereo(int16_t* l, int16_t* r, int count) {
    for (int i = 0; i < count; i++) {
        int32_t M = l[i], S = r[i];
        l[i] = clamp_q15(M + S);
        r[i] = clamp_q15(M - S);
    }
}

// Encoder
int encodeMADPCMBlock(const int16_t* inSamples, uint8_t* outBuffer, int16_t* history) {
    MADPCMBlockHeader* header = (MADPCMBlockHeader*)outBuffer;
    uint8_t* payload = outBuffer + sizeof(MADPCMBlockHeader);

    // High Precision Analysis
    int64_t r[LPC_ORDER + 1] = { 0 };
    computeAutocorr(inSamples, BLOCK_SAMPLES, LPC_ORDER, r);

    if (!computeLPC(r, LPC_ORDER, header->coeffs)) {
        memset(header->coeffs, 0, sizeof(header->coeffs));
    }

    // Simulation & Heuristics Loop (Preserved)
    int64_t sigEnergy = 0, resEnergy = 0;
    int zeroCrossings = 0;
    int32_t maxResSub[4] = { 0 };
    int32_t maxSigSub[4] = { 0 };
    int16_t simHist[LPC_ORDER];
    memcpy(simHist, history, sizeof(simHist));

    int32_t c0 = header->coeffs[0], c1 = header->coeffs[1];
    int32_t c2 = header->coeffs[2], c3 = header->coeffs[3];
    int badModel = 0;

    for (int i = 0; i < BLOCK_SAMPLES; i++) {
        int subIdx = i >> 8;

        // ZCR Check: (current ^ prev) < 0 means sign bit changed
        if (i > 0 && ((inSamples[i] ^ inSamples[i - 1]) < 0))
            zeroCrossings++;

        int32_t pred = (c0 * simHist[0]) + (c1 * simHist[1]) + (c2 * simHist[2]) + (c3 * simHist[3]);
        pred >>= 12;

        int32_t res = inSamples[i] - pred;
        sigEnergy += (int64_t)inSamples[i] * inSamples[i];
        resEnergy += (int64_t)res * res;

        if (abs_int32(res) > maxResSub[subIdx]) maxResSub[subIdx] = abs_int32(res);
        if (abs_int32(inSamples[i]) > maxSigSub[subIdx]) maxSigSub[subIdx] = abs_int32(inSamples[i]);

        if (abs_int32(res) > 28000) badModel = 1;

        simHist[3] = simHist[2]; simHist[2] = simHist[1]; simHist[1] = simHist[0]; simHist[0] = inSamples[i];
    }

    int highZCR = zeroCrossings > (BLOCK_SAMPLES >> 2);
    int lowGain = sigEnergy < (8 * (resEnergy + 1));
    int useFallback = badModel || (highZCR && lowGain);

    if (useFallback) {
        memset(header->coeffs, 0, sizeof(header->coeffs));
        c0 = c1 = c2 = c3 = 0;
    }

    // Populate scales with Gaussian Correction (~0.66x)
    for (int k = 0; k < 4; k++) {
        int32_t targetMax = useFallback ? maxSigSub[k] : maxResSub[k];
        // Gauss table max is ~1.5x, so we scale down to fit the peak
        int32_t s = (targetMax * 2) / 3;
        if (s < 10) s = 10;
        header->scales[k] = clamp_q15(s);
    }

    // Trellis Quantization (Replaces Closed Loop)
    int16_t curHist[LPC_ORDER];
    memcpy(curHist, history, sizeof(curHist));

    // Process 4 sub-blocks
    for (int sb = 0; sb < 4; sb++) {
        int16_t scale = header->scales[sb];
        const int16_t* blk = inSamples + (sb * 256);

        TrellisPath paths[TRELLIS_K];
        int activePaths = 1;
        paths[0].err = 0;
        memcpy(paths[0].h, curHist, sizeof(curHist));

        for (int i = 0; i < 256; i++) {
            TrellisCandidate cands[TRELLIS_K * 3];
            int candCount = 0;

            for (int p = 0; p < activePaths; p++) {
                int16_t* h = paths[p].h;
                int32_t pred = (c0 * h[0] + c1 * h[1] + c2 * h[2] + c3 * h[3]) >> 12;
                int32_t target = blk[i] - pred;

                // Check best greedy + neighbors
                int8_t bestN = findBestNibble(target, scale);
                int8_t centerIdx = (bestN & 0x08) ? (bestN & 0x07) + 8 : bestN;

                for (int offset = -1; offset <= 1; offset++) {
                    int idx = centerIdx + offset;
                    if (idx < 0 || idx > 15) continue;

                    int8_t nib = (idx < 8) ? idx : (idx & 0x07) | 0x08;
                    int16_t rec = clamp_q15(pred + dequantizeLUT(nib, scale));
                    int64_t sqErr = (int64_t)(blk[i] - rec) * (blk[i] - rec);

                    cands[candCount].err = paths[p].err + sqErr;
                    cands[candCount].pIdx = p;
                    cands[candCount].nib = nib;
                    cands[candCount].val = rec;
                    candCount++;
                }
            }

            // Select best K
            TrellisPath nextPaths[TRELLIS_K];
            int nextCount = 0;
            for (int k = 0; k < TRELLIS_K && k < candCount; k++) {
                int bestJ = -1;
                int64_t minE = -1;
                for (int j = 0; j < candCount; j++) {
                    if (cands[j].err == -1) continue;
                    if (bestJ == -1 || cands[j].err < minE) { minE = cands[j].err; bestJ = j; }
                }
                if (bestJ == -1) break;

                nextPaths[k] = paths[cands[bestJ].pIdx];
                nextPaths[k].err = cands[bestJ].err;
                nextPaths[k].path[i] = cands[bestJ].nib;

                int16_t* h = nextPaths[k].h;
                h[3] = h[2]; h[2] = h[1]; h[1] = h[0]; h[0] = cands[bestJ].val;

                cands[bestJ].err = -1; // Mark used
                nextCount++;
            }
            activePaths = nextCount;
            memcpy(paths, nextPaths, sizeof(TrellisPath) * activePaths);
        }

        // Commit best path
        memcpy(curHist, paths[0].h, sizeof(curHist));
        int byteStart = sb * 128;
        for (int i = 0; i < 256; i += 2) {
            payload[byteStart + (i >> 1)] = (paths[0].path[i] & 0x0F) | (paths[0].path[i + 1] << 4);
        }
    }

    memcpy(history, curHist, sizeof(curHist));
    return sizeof(MADPCMBlockHeader) + (BLOCK_SAMPLES / 2);
}

// Adaptive Stereo Encoder
int encodeMADPCMStereoBlock(const int16_t* inL, const int16_t* inR, uint8_t* outBuffer, int16_t* hL, int16_t* hR) {
    int16_t bufL[BLOCK_SAMPLES], bufR[BLOCK_SAMPLES];
    memcpy(bufL, inL, sizeof(bufL));
    memcpy(bufR, inR, sizeof(bufR));

    // Adaptive Check: Calculate Energies
    int64_t eLR = 0, eMS = 0;
    for (int i = 0; i < BLOCK_SAMPLES; i += 4) {
        int32_t m = (bufL[i] + bufR[i]) >> 1;
        int32_t s = (bufL[i] - bufR[i]) >> 1;
        eLR += abs_int32(bufL[i]) + abs_int32(bufR[i]);
        eMS += abs_int32(m) + abs_int32(s);
    }

    bool useMS = (eMS < eLR * 0.9f);
    if (useMS) stereoToMidSide(bufL, bufR, BLOCK_SAMPLES);

    int s1 = encodeMADPCMBlock(bufL, outBuffer, hL);
    int s2 = encodeMADPCMBlock(bufR, outBuffer + s1, hR);

    // Flagging: Use LSB of first scale
    MADPCMBlockHeader* h = (MADPCMBlockHeader*)outBuffer;
    if (useMS) h->scales[0] |= 1;
    else       h->scales[0] &= ~1;

    return s1 + s2;
}

// Decoder
void decodeMADPCMBlock(const uint8_t* inBuffer, int16_t* outSamples, int16_t* history)
{
    const MADPCMBlockHeader* header = (const MADPCMBlockHeader*)inBuffer;
    const uint8_t* payload = inBuffer + sizeof(MADPCMBlockHeader);

    const int32_t c0 = header->coeffs[0];
    const int32_t c1 = header->coeffs[1];
    const int32_t c2 = header->coeffs[2];
    const int32_t c3 = header->coeffs[3];

    int16_t h0 = history[0];
    int16_t h1 = history[1];
    int16_t h2 = history[2];
    int16_t h3 = history[3];

    int sample = 0;

    // Each scale applies to 256 samples
    for (int s = 0; s < BLOCK_SAMPLES / 256; s++) {
        int16_t scale = header->scales[s] & ~1;

        for (int i = 0; i < 256; i++, sample++) {
            int32_t pred =
                (c0 * h0 + c1 * h1 + c2 * h2 + c3 * h3) >> 12;

            uint8_t byte = payload[sample >> 1];
            int8_t nibble = (sample & 1)
                ? (byte >> 4)
                : (byte & 0x0F);

            int16_t delta = dequantizeLUT(nibble, scale);
            int32_t val = pred + delta;

            h3 = h2;
            h2 = h1;
            h1 = h0;
            h0 = clamp_q15(val);

            outSamples[sample] = h0;
        }
    }

    history[0] = h0; history[1] = h1; history[2] = h2; history[3] = h3;
}

void decodeMADPCMStereoBlock(const uint8_t* inBuffer, int16_t* outL, int16_t* outR, int16_t* hL, int16_t* hR) {
    const MADPCMBlockHeader* h = (const MADPCMBlockHeader*)inBuffer;
    bool isMS = (h->scales[0] & 1);

    int sz = sizeof(MADPCMBlockHeader) + (BLOCK_SAMPLES / 2);
    decodeMADPCMBlock(inBuffer, outL, hL);
    decodeMADPCMBlock(inBuffer + sz, outR, hR);

    if (isMS) midSideToStereo(outL, outR, BLOCK_SAMPLES);
}