# Modern Adaptive Differential Pulse-Code Modulation

ADPCM didn't get worse. Everyone just stopped trying.  
A small, fast, waveform audio codec that fixes what IMA ADPCM got wrong.  
No transforms. No psychoacoustics. No entropy coding. No excuses.

“The obvious replacement for IMA ADPCM when you control both ends”

## what is it?

The math of FLAC, but the constraints of ADPCM.  
~500 lines of pure C.

## under the hood

MADPCM implements a fixed-point 4th-order Linear Predictive Coder (LPC) derived via Levinson-Durbin recursion.  
A pre-encoding simulation and heuristics loop validates coefficients using Zero Crossing Rate (ZCR) and signal-to-residual energy ratios, triggering fallback mechanisms on model instability.  
Quantization utilizes a Viterbi trellis search ($K=4$) for global optimization or a greedy O(1) inverse lookahead for speed.  
Stereo coupling employs adaptive Mid-Side coding, dynamically switching based on inter-channel energy correlation.  
The core uses Gaussian-optimized step tables and is strictly integer-based to support toasters (read: freestanding environments).

The decoder runs entirely in registers, no stack spilling provided it's compiled with optimizations (tested with mingw gcc 15.1.0).

## warning

This is work in progress, expect breaking changes every other commit.

## the numbers

Measured on 16-bit PCM sources:

**IMA ADPCM** (ffmpeg):  
| File                       | RMSE  | MAE  | PSNR     | SNR      | Enc Time | Dec Time |
| -------------------------- | ----- | ---- | -------- | -------- | -------- | -------- |
| Marche_Persanne.wav        | 49.81 | 3345 | 56.36 dB | 32.38 dB | 0.119s   | 0.058s   |
| Egyptischer_Marsch.wav     | 56.51 | 5872 | 55.27 dB | 32.70 dB | 0.256s   | 0.145s   |
| suppe_poet_and_peasant.wav | 65.80 | 8322 | 53.94 dB | 31.26 dB | 0.636s   | 0.348s   |
| the_four_seasons.wav       | 47.71 | 2630 | 56.74 dB | 36.19 dB | 3.927s   | 2.150s   |

**ADPCM-XQ** (lookahead=8):  
| File                       | RMSE  | MAE  | PSNR     | SNR      | Enc Time |
| -------------------------- | ----- | ---- | -------- | -------- | -------- |
| Marche_Persanne.wav        | 42.75 | 1228 | 57.69 dB | 33.71 dB | 17.486s  |
| Egyptischer_Marsch.wav     | 47.31 | 1851 | 56.81 dB | 34.25 dB | 40.850s  |
| suppe_poet_and_peasant.wav | 54.54 | 1640 | 55.58 dB | 32.89 dB | 102.432s |
| the_four_seasons.wav       | 44.67 | 778  | 57.31 dB | 36.76 dB | 607.841s |

**MADPCM** (slow):  
| File                       | RMSE  | MAE  | PSNR     | SNR      | Enc Time | Dec Time |
| -------------------------- | ----- | ---- | -------- | -------- | -------- | -------- |
| Marche_Persanne.wav        | 27.30 | 343  | 61.59 dB | 37.61 dB | 0.178s   | 0.008s   |
| Egyptischer_Marsch.wav     | 32.09 | 623  | 60.18 dB | 37.62 dB | 0.384s   | 0.016s   |
| suppe_poet_and_peasant.wav | 37.05 | 810  | 58.93 dB | 36.25 dB | 0.991s   | 0.038s   |
| the_four_seasons.wav       | 25.00 | 222  | 62.35 dB | 41.80 dB | 6.065s   | 0.214s   |

**MADPCM** (fast; no lookahead):  
| File                       | RMSE  | MAE  | PSNR     | SNR      | Enc Time |
| -------------------------- | ----- | ---- | -------- | -------- | -------- |
| Marche_Persanne.wav        | 29.01 | 512  | 61.06 dB | 37.08 dB | 0.020s   |
| Egyptischer_Marsch.wav     | 34.01 | 600  | 59.68 dB | 37.11 dB | 0.054s   |
| suppe_poet_and_peasant.wav | 39.36 | 802  | 58.41 dB | 35.72 dB | 0.115s   |
| the_four_seasons.wav       | 26.37 | 271  | 61.89 dB | 41.34 dB | 0.680s   |

All tested on a Ryzen 5 3600XT running Windows Server 2025 Datacenter build 26100.32230  
Decode speed for MADPCM is the same on slow and fast, this affects only the encoder.  
Encode speed for ADPCM-XQ includes file I/O, tested on a ramdisk to minimize potential penalty.  
File size reduction: ~74% vs WAV, 4.48% overhead compared to IMA ADPCM.

<details><summary>test file details</summary>

All downloaded from youtube with premium at max quality, converted to wav.  
Egyptischer_Marsch.wav 43,5 MB Length: 03:57 Sample rate: 48 kHz Sample size: 16 bit Stereo  
Marche_Persanne.wav 18,7 MB Length: 01:51 Sample rate: 48 kHz Sample size: 16 bit Stereo  
suppe_poet_and_peasant.wav 107 MB Length: 09:49 Sample rate: 48 kHz Sample size: 16 bit Stereo  
the_four_seasons.wav 665 MB Length: 01:00:33 Sample rate: 48 kHz Sample size: 16 bit Stereo

</details>

<details><summary>EBU Sound Quality Assessment Material</summary>

**MADPCM** (slow):  
| Audio Sample  | RMSE    | MAE  | PSNR  | SNR   |
| ------------- | ------- | ---- | ----- | ----- |
| Trumpet       | 16.0853 | 259  | 66.18 | 41.73 |
| Castanets     | 80.8316 | 2502 | 52.16 | 20.94 |
| Side Drum     | 16.3826 | 317  | 66.02 | 39.33 |
| Glockenspiel  | 45.8915 | 1782 | 57.07 | 30.78 |
| Vibraphone    | 17.9936 | 256  | 65.21 | 42.56 |
| Grand Piano   | 15.4270 | 213  | 66.54 | 43.23 |
| Harpsichord   | 30.6090 | 770  | 60.59 | 27.78 |
| Female Speech | 45.9353 | 1135 | 57.07 | 33.53 |

**IMA ADPCM** (ffmpeg):  
| Audio Sample  | RMSE    | MAE   | PSNR  | SNR   |
| ------------- | ------- | ----- | ----- | ----- |
| Trumpet       | 34.8689 | 2434  | 59.46 | 35.01 |
| Castanets     | 109.273 | 15673 | 49.54 | 18.33 |
| Side Drum     | 24.5968 | 2407  | 62.49 | 35.80 |
| Glockenspiel  | 83.6360 | 20459 | 51.86 | 25.57 |
| Vibraphone    | 23.4344 | 865   | 62.91 | 40.26 |
| Grand Piano   | 16.6159 | 1233  | 65.90 | 42.58 |
| Harpsichord   | 39.5648 | 3725  | 58.36 | 25.55 |
| Female Speech | 72.1851 | 3461  | 53.14 | 29.61 |

</details>

<details><summary>how other underpowered hardware behaves</summary>

samsung s21:  
```
~/madpcm $ ./test.out encode slow suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
Loaded suppe_poet_and_peasant.wav: 48000 Hz, Stereo
Encoded in 4.128757 seconds (slow mode)
~/madpcm $ ./test.out encode fast suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
Loaded suppe_poet_and_peasant.wav: 48000 Hz, Stereo
Encoded in 0.311960 seconds (fast mode)
~/madpcm $ ./test.out decode suppe_poet_and_peasant_ad.wav suppe_poet_and_peasant_d.wav
Loaded suppe_poet_and_peasant_ad.wav: 2 channels, 48000 Hz, 28272640 samples
Decoded in 0.079933 seconds
```

raspberry pi 4:  
```
maksw@raspberrypi:~/madpcm $ ./test.out encode slow suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
Loaded suppe_poet_and_peasant.wav: 48000 Hz, Stereo
Encoded in 4.307982 seconds (slow mode)
maksw@raspberrypi:~/madpcm $ ./test.out encode fast suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
Loaded suppe_poet_and_peasant.wav: 48000 Hz, Stereo
Encoded in 0.430510 seconds (fast mode)
maksw@raspberrypi:~/madpcm $ ./test.out decode suppe_poet_and_peasant_ad.wav suppe_poet_and_peasant_d.wav
Loaded suppe_poet_and_peasant_ad.wav: 2 channels, 48000 Hz, 28272640 samples
Decoded in 0.213578 seconds
```

</details>

<details><summary>how a 486 behaves (dosbox-x)</summary>

**Environment:** Emulated 486 (25,000 cycles), 15s Stereo WAV @ 44.1kHz.

```text
C:\>mad.exe encode fast 37.wav 37_ad.wav
Loaded 37.wav: 44100 Hz, Stereo
Encoded in 10.659341 seconds (fast mode)

C:\>mad.exe decode 37_ad.wav 37_d.wav
Loaded 37_ad.wav: 2 channels, 44100 Hz, 661504 samples
Decoded in 2.747253 seconds
```

Encoder runs faster than realtime.  
Decoder is ~5.4× realtime.  
scary shit

</details>

## does it sound better?

To me. I think. Better than IMA ADPCM.  
No formal ABX tests were done yet.

In my personal testing (Sony WF-1000XM4 connected via LDAC) compared to the source, MADPCM adds a very small, mostly uncorrelated noise floor. Compared to IMA ADPCM, it sounds much cleaner.

## how 2 use it

it's an `stb` style library:  
include `madpcm.h` in your code.  
`MADPCM_IMPLEMENTATION` for the actual implementation (define **once**)  
`MADPCM_FREESTANDING` builds the engine for freestanding  
`MADPCM_MEMFUNCS` defines own memcpy and memset  
Use when MSVC /Oi- or always define it for freestanding gcc/clang ***if*** you don't have these functions.

## what this does not care about

- standards
- legacy compatibility
- your favorite codec comparison

## how to contribute

Test it. Report bugs. Open PRs.  
4 spaces, K&R styling.

## q&a

Q: Why make this?  
A: IMA ADPCM works. This works better. IMA's everywhere because it's cheap and simple, not because it's good.

Q: Is this better than Opus/AAC/MP3?  
A: If you squint hard enough and just look at the numbers. This is a scalpel, not a swiss army knife.

Q: Is this a competitor to Opus?  
A: No.

Q: What's the use case?  
A: Smaller-than-WAV files, predictable latency, trivial decoding cost.

Q: Are there transforms or psychoacoustic models?  
A: No. Conceptually inspired by lossless predictive codecs like FLAC.

Q: Why is encoding slow?  
A: Quality costs time. Decode is what matters for playback.

Q: Will you add [feature]?  
A: Depends. Open an issue.

Q: Can I use this in production?  
A: No warranty. Assume breaking changes until it's finished.  
If you *really* want to use it, at least keep track of the commit you used.
