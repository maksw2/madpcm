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
static inline int16_t madpcm__clamp_q15(int32_t value) {
    if (value > MADPCM_Q15_MAX) return MADPCM_Q15_MAX;
    if (value < MADPCM_Q15_MIN) return MADPCM_Q15_MIN;
    return (int16_t)value;
}
static inline int32_t madpcm__abs_int32(int32_t value) { return (value < 0) ? -value : value; }
static inline int64_t madpcm__abs_int64(int64_t value) { return (value < 0) ? -value : value; }

static void madpcm__compute_autocorr(const int16_t* signal, int length, int order, int64_t* r_out) {
    for (int i = 0; i <= order; i++) {
        int64_t sum = 0;
        for (int j = 0; j < length - i; j++) sum += (int64_t)signal[j] * signal[j + i];
        r_out[i] = sum;
    }
    r_out[0] = (r_out[0] > 256) ? r_out[0] : 256;
    r_out[0] += (r_out[0] >> 9);
}

static bool madpcm__compute_LPC(const int64_t* r, int order, int16_t* outCoeffs) {
    int32_t a[MADPCM_LPC_ORDER + 1][MADPCM_LPC_ORDER + 1] = { 0 };
    int64_t error[MADPCM_LPC_ORDER + 1];
    error[0] = r[0];

    for (int k = 1; k <= order; k++) {
        int64_t sum = 0;
        for (int j = 1; j < k; j++) {
            int64_t product = (int64_t)a[k - 1][j] * r[k - j];
            sum += (product >> 24);
        }
        int64_t numerator = r[k] - sum;
        int64_t denominator = error[k - 1];
        if (denominator == 0) denominator = 1;
        if (madpcm__abs_int64(numerator) >= madpcm__abs_int64(denominator)) return false;

        int32_t reflection_coeff = (int32_t)((numerator << 24) / denominator);
        a[k][k] = reflection_coeff;

        for (int j = 1; j < k; j++) {
            int64_t term = (int64_t)reflection_coeff * a[k - 1][k - j];
            a[k][j] = a[k - 1][j] - (int32_t)(term >> 24);
        }

        int64_t reflection_sq = (int64_t)reflection_coeff * reflection_coeff;
        error[k] = (error[k - 1] * ((1LL << 30) - (reflection_sq >> 18))) >> 30;
    }

    for (int i = 0; i < order; i++) {
        outCoeffs[i] = madpcm__clamp_q15(a[order][i + 1] >> 12);
    }
    return true;
}

// LUT-based Dequantizer
static inline int16_t madpcm__dequantize_LUT(int8_t nibble, int16_t scale) {
    int8_t index = nibble & 0x07;
    int8_t sign = (nibble & 0x08) ? -1 : 1;
    int32_t value = ((int32_t)scale * madpcm__gauss_steps[index]) >> 11;
    return (sign < 0) ? (int16_t)-value : (int16_t)value;
}

// Finds best nibble in O(1) using inverse quantization estimation
static inline int8_t madpcm__find_best_nibble(int32_t target, int16_t scale) {
    if (scale < 1) scale = 1;

    int32_t abs_target = madpcm__abs_int32(target);
    // Approximation: index ~ (val * 160) >> 16
    // 1/410 (step size) in Q16 is approx 160
    int32_t normalized_val = (abs_target << 11) / scale;
    int32_t index = (normalized_val * 160) >> 16;

    if (index > 7) index = 7;

    // Check neighbor to handle boundary cases
    int32_t estimate = (scale * madpcm__gauss_steps[index]) >> 11;
    int32_t diff1 = madpcm__abs_int32(abs_target - estimate);

    if (index < 7) {
        int32_t estimate_next = (scale * madpcm__gauss_steps[index + 1]) >> 11;
        int32_t diff2 = madpcm__abs_int32(abs_target - estimate_next);
        if (diff2 < diff1) index++;
    }

    if (target < 0) return (int8_t)(index | 0x08);
    return (int8_t)index;
}

static inline void madpcm__stereo_to_mid_side(int16_t* left, int16_t* right, int count) {
    for (int i = 0; i < count; i++) {
        int32_t l_val = left[i], r_val = right[i];
        left[i] = (int16_t)((l_val + r_val) >> 1);
        right[i] = (int16_t)((l_val - r_val) >> 1);
    }
}

static inline void madpcm__mid_side_to_stereo(int16_t* left, int16_t* right, int count) {
    for (int i = 0; i < count; i++) {
        int32_t mid = left[i], side = right[i];
        left[i] = madpcm__clamp_q15(mid + side);
        right[i] = madpcm__clamp_q15(mid - side);
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
    int64_t autocorrelation[MADPCM_LPC_ORDER + 1] = { 0 };
    madpcm__compute_autocorr(inSamples, MADPCM_BLOCK_SAMPLES, MADPCM_LPC_ORDER, autocorrelation);

    if (!madpcm__compute_LPC(autocorrelation, MADPCM_LPC_ORDER, header->coeffs)) {
        madpcm_memset(header->coeffs, 0, sizeof(header->coeffs));
    }

    // Simulation & Heuristics Loop
    int64_t signal_energy = 0, residual_energy = 0;
    int zero_crossings = 0;
    int32_t max_residual_sub[4] = { 0 };
    int32_t max_signal_sub[4] = { 0 };
    
    // Initialize simulation history with perfect seeds
    int16_t sim_history[MADPCM_LPC_ORDER];
    madpcm_memcpy(sim_history, header->seed, sizeof(sim_history));

    int32_t c0 = header->coeffs[0], c1 = header->coeffs[1];
    int32_t c2 = header->coeffs[2], c3 = header->coeffs[3];
    bool bad_model = false;

    for (int i = 0; i < MADPCM_BLOCK_SAMPLES; i++) {
        int sub_block_idx = i >> 8;

        // ZCR Check: (current ^ prev) < 0 means sign bit changed
        if (i > 0 && ((inSamples[i] ^ inSamples[i - 1]) < 0))
            zero_crossings++;

        int32_t prediction = (c0 * sim_history[0]) + (c1 * sim_history[1]) + (c2 * sim_history[2]) + (c3 * sim_history[3]);
        prediction >>= 12;

        int32_t residual = inSamples[i] - prediction;
        signal_energy += (int64_t)inSamples[i] * inSamples[i];
        residual_energy += (int64_t)residual * residual;

        if (madpcm__abs_int32(residual) > max_residual_sub[sub_block_idx]) max_residual_sub[sub_block_idx] = madpcm__abs_int32(residual);
        if (madpcm__abs_int32(inSamples[i]) > max_signal_sub[sub_block_idx]) max_signal_sub[sub_block_idx] = madpcm__abs_int32(inSamples[i]);

        if (madpcm__abs_int32(residual) > 28000) bad_model = true;

        sim_history[3] = sim_history[2]; sim_history[2] = sim_history[1]; sim_history[1] = sim_history[0]; sim_history[0] = inSamples[i];
    }

    bool high_zcr = zero_crossings > (MADPCM_BLOCK_SAMPLES >> 2);
    bool low_gain = signal_energy < (8 * (residual_energy + 1));
    bool use_fallback = bad_model || (high_zcr && low_gain);

    if (use_fallback) {
        madpcm_memset(header->coeffs, 0, sizeof(header->coeffs));
        c0 = c1 = c2 = c3 = 0;
    }

    // Populate scales with Gaussian Correction (~0.66x)
    for (int k = 0; k < 4; k++) {
        int32_t target_max = use_fallback ? max_signal_sub[k] : max_residual_sub[k];
        int32_t scale_val = (target_max * 2) / 3;
        if (scale_val < 10) scale_val = 10;
        header->scales[k] = madpcm__clamp_q15(scale_val);
    }

    // Trellis Quantization
    int16_t current_history[MADPCM_LPC_ORDER];
    madpcm_memcpy(current_history, header->seed, sizeof(current_history));

    // Process 4 sub-blocks
    for (int sb = 0; sb < 4; sb++) {
        int16_t current_scale = header->scales[sb];
        const int16_t* block_ptr = inSamples + (sb * 256);
        uint8_t* sub_payload_ptr = payload + (sb * 128);

        // FAST MODE: No Trellis
        if (fast) {
            for (int i = 0; i < 256; i += 2) {
                // Nibble 1
                int32_t p1 = (c0 * current_history[0] + c1 * current_history[1] + c2 * current_history[2] + c3 * current_history[3]) >> 12;
                int8_t n1 = madpcm__find_best_nibble(block_ptr[i] - p1, current_scale);
                int16_t r1 = madpcm__clamp_q15(p1 + madpcm__dequantize_LUT(n1, current_scale));
                current_history[3] = current_history[2]; current_history[2] = current_history[1]; current_history[1] = current_history[0]; current_history[0] = r1;

                // Nibble 2
                int32_t p2 = (c0 * current_history[0] + c1 * current_history[1] + c2 * current_history[2] + c3 * current_history[3]) >> 12;
                int8_t n2 = madpcm__find_best_nibble(block_ptr[i + 1] - p2, current_scale);
                int16_t r2 = madpcm__clamp_q15(p2 + madpcm__dequantize_LUT(n2, current_scale));
                current_history[3] = current_history[2]; current_history[2] = current_history[1]; current_history[1] = current_history[0]; current_history[0] = r2;

                sub_payload_ptr[i >> 1] = (n1 & 0x0F) | (n2 << 4);
            }
            continue;
        }

        // QUALITY MODE: Trellis
        // Traceback buffers [sample][path]
        int8_t trace_nibble[256][MADPCM_TRELLIS_K];
        int8_t trace_index[256][MADPCM_TRELLIS_K];

        madpcm__trellis_state current_states[MADPCM_TRELLIS_K];
        int active_paths = 1;
        current_states[0].err = 0;
        madpcm_memcpy(current_states[0].h, current_history, sizeof(current_history));

        for (int sample_idx = 0; sample_idx < 256; sample_idx++) {
            madpcm__trellis_candidate candidates[MADPCM_TRELLIS_K * 3];
            int candidate_count = 0;

            for (int path_idx = 0; path_idx < active_paths; path_idx++) {
                int16_t* path_hist = current_states[path_idx].h;
                int32_t prediction = (c0 * path_hist[0] + c1 * path_hist[1] + c2 * path_hist[2] + c3 * path_hist[3]) >> 12;
                int32_t residual_target = block_ptr[sample_idx] - prediction;

                // Check best greedy + neighbors
                int8_t best_nibble = madpcm__find_best_nibble(residual_target, current_scale);
                int8_t center_idx = (best_nibble & 0x08) ? (best_nibble & 0x07) + 8 : best_nibble;

                for (int offset = -1; offset <= 1; offset++) {
                    int search_idx = center_idx + offset;
                    if (search_idx < 0 || search_idx > 15) continue;

                    int8_t current_nib = (search_idx < 8) ? (int8_t)search_idx : (int8_t)((search_idx & 0x07) | 0x08);
                    int16_t reconstructed = madpcm__clamp_q15(prediction + madpcm__dequantize_LUT(current_nib, current_scale));
                    int64_t diff = block_ptr[sample_idx] - reconstructed;
                    int64_t squared_error = diff * diff;

                    candidates[candidate_count].err = current_states[path_idx].err + squared_error;
                    candidates[candidate_count].p_idx = path_idx;
                    candidates[candidate_count].nib = current_nib;
                    candidates[candidate_count].val = reconstructed;
                    candidate_count++;
                }
            }

            // Select best K
            madpcm__trellis_state next_states[MADPCM_TRELLIS_K];
            int next_count = 0;

            for (int k = 0; k < MADPCM_TRELLIS_K && k < candidate_count; k++) {
                int best_cand_idx = -1;
                int64_t min_err = -1;
                for (int j = 0; j < candidate_count; j++) {
                    if (candidates[j].err == -1) continue;
                    if (best_cand_idx == -1 || candidates[j].err < min_err) { min_err = candidates[j].err; best_cand_idx = j; }
                }
                if (best_cand_idx == -1) break;

                // Save traceback data instead of full path copy
                trace_index[sample_idx][k] = (int8_t)candidates[best_cand_idx].p_idx;
                trace_nibble[sample_idx][k] = candidates[best_cand_idx].nib;

                next_states[k] = current_states[candidates[best_cand_idx].p_idx];
                next_states[k].err = candidates[best_cand_idx].err;

                int16_t* next_hist = next_states[k].h;
                next_hist[3] = next_hist[2]; next_hist[2] = next_hist[1]; next_hist[1] = next_hist[0]; next_hist[0] = candidates[best_cand_idx].val;

                candidates[best_cand_idx].err = -1; // Mark used
                next_count++;
            }
            active_paths = next_count;
            madpcm_memcpy(current_states, next_states, sizeof(madpcm__trellis_state) * active_paths);
        }

        // Commit best path (Backtracking)
        madpcm_memcpy(current_history, current_states[0].h, sizeof(current_history));
        int backtrack_path_idx = 0; // Index 0 is best due to sort

        uint8_t final_nibbles[256];
        for (int i = 255; i >= 0; i--) {
            final_nibbles[i] = trace_nibble[i][backtrack_path_idx];
            backtrack_path_idx = trace_index[i][backtrack_path_idx];
        }

        for (int j = 0; j < 128; j++) {
            sub_payload_ptr[j] = (final_nibbles[2 * j] & 0x0F) | (final_nibbles[2 * j + 1] << 4);
        }
    }

    return sizeof(madpcm__block_header) + (MADPCM_BLOCK_SAMPLES / 2);
}

// Adaptive Stereo Encoder
int madpcm_encode_block_stereo(const int16_t* inL, const int16_t* inR, const int16_t* prevL, const int16_t* prevR, uint8_t* outBuffer, bool fast) {
    int16_t buffer_left[MADPCM_BLOCK_SAMPLES], buffer_right[MADPCM_BLOCK_SAMPLES];
    int16_t history_left[4], history_right[4];

    madpcm_memcpy(buffer_left, inL, sizeof(buffer_left));
    madpcm_memcpy(buffer_right, inR, sizeof(buffer_right));

    // Prepare local history copies (safe to modify)
    if (prevL && prevR) {
        madpcm_memcpy(history_left, prevL, 4 * sizeof(int16_t));
        madpcm_memcpy(history_right, prevR, 4 * sizeof(int16_t));
    } else {
        madpcm_memset(history_left, 0, sizeof(history_left));
        madpcm_memset(history_right, 0, sizeof(history_right));
    }

    // Adaptive Check: Calculate Energies
    int64_t energy_lr = 0, energy_ms = 0;
    for (int i = 0; i < MADPCM_BLOCK_SAMPLES; i += 4) {
        int32_t mid = (buffer_left[i] + buffer_right[i]) >> 1;
        int32_t side = (buffer_left[i] - buffer_right[i]) >> 1;
        energy_lr += madpcm__abs_int32(buffer_left[i]) + madpcm__abs_int32(buffer_right[i]);
        energy_ms += madpcm__abs_int32(mid) + madpcm__abs_int32(side);
    }

    bool use_ms = (energy_ms * 10 < energy_lr * 9);
    
    // If using MS, we must also convert the HISTORY to MS
    if (use_ms) {
        madpcm__stereo_to_mid_side(buffer_left, buffer_right, MADPCM_BLOCK_SAMPLES);
        madpcm__stereo_to_mid_side(history_left, history_right, 4); 
    }

    int size_left = madpcm_encode_block(buffer_left, history_left, outBuffer, fast);
    int size_right = madpcm_encode_block(buffer_right, history_right, outBuffer + size_left, fast);

    // Flagging: Use LSB of first scale
    madpcm__block_header* main_header = (madpcm__block_header*)outBuffer;
    main_header->scales[0] = (main_header->scales[0] & ~1) | (use_ms ? 1 : 0);

    return size_left + size_right;
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

    int current_sample_idx = 0;

    // Each scale applies to 256 samples
    for (int s_idx = 0; s_idx < MADPCM_BLOCK_SAMPLES / 256; s_idx++) {
        int16_t current_scale = header->scales[s_idx] & ~1;

        for (int i = 0; i < 256; i++, current_sample_idx++) {
            int32_t prediction = (c0 * h0 + c1 * h1 + c2 * h2 + c3 * h3) >> 12;

            uint8_t payload_byte = payload[current_sample_idx >> 1];
            int8_t current_nibble = (current_sample_idx & 1)
                ? (payload_byte >> 4)
                : (payload_byte & 0x0F);

            int16_t delta = madpcm__dequantize_LUT(current_nibble, current_scale);
            int32_t final_val = prediction + delta;

            h3 = h2;
            h2 = h1;
            h1 = h0;
            h0 = madpcm__clamp_q15(final_val);

            outSamples[current_sample_idx] = h0;
        }
    }
}

void madpcm_decode_block_stereo(const uint8_t* inBuffer, int16_t* outL, int16_t* outR) {
    const madpcm__block_header* header = (const madpcm__block_header*)inBuffer;
    bool is_mid_side = (header->scales[0] & 1);

    int mono_block_size = sizeof(madpcm__block_header) + (MADPCM_BLOCK_SAMPLES / 2);
    madpcm_decode_block(inBuffer, outL);
    madpcm_decode_block(inBuffer + mono_block_size, outR);

    if (is_mid_side) madpcm__mid_side_to_stereo(outL, outR, MADPCM_BLOCK_SAMPLES);
}

#endif // MADPCM_IMPLEMENTATION
#endif // MADPCM_H