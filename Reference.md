# MADPCM Implementation Reference

### Compile-time Defines

| Define | Effect |
|---|---|
| `MADPCM_IMPLEMENTATION` | Emits the implementation. Define in one TU only. |
| `MADPCM_FREESTANDING` | Strips `<string.h>` dependency. You must supply memcpy/memset. |
| `MADPCM_INTERNAL_MEMFUNCS` | Uses built-in memcpy/memset implementations. Use with `MADPCM_FREESTANDING` when the runtime doesn't provide them. |
| `MADPCM_MEMSET` / `MADPCM_MEMCPY` | Bring-your-own replacements. Both must be defined together. Only valid under `MADPCM_FREESTANDING`. |

---

## Constants

```c
#define MADPCM_LPC_ORDER      4      // predictor taps
#define MADPCM_BLOCK_SAMPLES  1024   // PCM samples per block (both channels)
#define MADPCM_HEADER_SIZE    24     // bytes — one mono block header
#define MADPCM_PAYLOAD_SIZE   512    // bytes — one mono block payload
#define MADPCM_MONO_SIZE      536    // HEADER + PAYLOAD
#define MADPCM_STEREO_SIZE    1072   // MONO_SIZE * 2
#define MADPCM_WAVE_FORMAT    0x4D41 // RIFF fmt chunk audio format tag
```

---

## API

### Encoding

```c
int madpcm_encode_block(
    const int16_t* inSamples,   // [MADPCM_BLOCK_SAMPLES] input PCM
    const int16_t* prevSamples, // [4] last 4 samples of previous block, or NULL
    uint8_t*       outBuffer,   // [MADPCM_MONO_SIZE] output
    bool           fast         // true = greedy O(1), false = Viterbi trellis
);

int madpcm_encode_block_stereo(
    const int16_t* inL,       // [MADPCM_BLOCK_SAMPLES]
    const int16_t* inR,       // [MADPCM_BLOCK_SAMPLES]
    const int16_t* prevL,     // [4] or NULL
    const int16_t* prevR,     // [4] or NULL
    uint8_t*       outBuffer, // [MADPCM_STEREO_SIZE]
    bool           fast
);
```

`prevSamples` / `prevL` / `prevR` seed the LPC history across block boundaries. Pass `NULL` for the first block; pass a pointer to the 4 samples immediately preceding the current block offset for all subsequent blocks (see `test.c` for the exact pattern). Skipping this will cause a small discontinuity at every block boundary.

Both functions return the number of bytes written, which will always equal `MADPCM_MONO_SIZE` or `MADPCM_STEREO_SIZE` respectively.

For the decoder: input buffer must be exactly MADPCM_MONO_SIZE or MADPCM_STEREO_SIZE bytes — there's no partial decode, you get 1024 samples or nothing. Output must fit MADPCM_BLOCK_SAMPLES int16_t values per channel.

For the encoder: input must be exactly MADPCM_BLOCK_SAMPLES samples. If your source audio doesn't divide evenly, you have to pad the last block to MADPCM_BLOCK_SAMPLES yourself before calling encode. Zero-padding is the obvious choice.

`fast = true` is the default choice for real-time or batch work. `fast = false` runs a Viterbi trellis search (K=4) and produces lower RMSE at roughly 5–10× the encoding cost.

### Decoding

```c
void madpcm_decode_block(
    const uint8_t* inBuffer,   // [MADPCM_MONO_SIZE]
    int16_t*       outSamples  // [MADPCM_BLOCK_SAMPLES]
);

void madpcm_decode_block_stereo(
    const uint8_t* inBuffer, // [MADPCM_STEREO_SIZE]
    int16_t*       outL,     // [MADPCM_BLOCK_SAMPLES]
    int16_t*       outR      // [MADPCM_BLOCK_SAMPLES]
);
```

Decode has no `fast` parameter as it's always the same path. Blocks are fully independent, you can decode them in parallel, in any order, from any offset. No cross-block state is needed on decode.

---

## Block Layout

A mono block is always `MADPCM_MONO_SIZE` (536) bytes.

```
[ madpcm__block_header : 24 bytes ][ payload : 512 bytes ]
```

## Parallelism

Blocks have no decode-time dependencies on each other. Encoding them in parallel is safe provided each thread writes to its own output region. `test.c` does this with `#pragma omp parallel for`; adapt as needed for your threading model.

---

## RIFF Container

If you're wrapping blocks in a WAV file:

| Field | Value |
|---|---|
| `wFormatTag` | `0x4D41` (`MADPCM_WAVE_FORMAT`) |
| `wBitsPerSample` | `0` (compressed) |
| `nBlockAlign` | `MADPCM_MONO_SIZE` or `MADPCM_STEREO_SIZE` |
| `nAvgBytesPerSec` | `(sampleRate / MADPCM_BLOCK_SAMPLES) * nBlockAlign` |

Standard players won't decode this. The format tag is unregistered. It's a container of convenience, not a compatibility claim.

---

## Porting Notes

- All math is fixed-point. No floats in the codec path.
- The codec itself allocates nothing.
- Stack usage in the slow encoder is non-trivial (trellis state arrays). On constrained targets, prefer `fast = true`.
- The decoder is designed to run entirely in registers with optimizations enabled. Verified with mingw gcc 15.1.0.
