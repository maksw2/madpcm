# Modern Adaptive Differential Pulse-Code Modulation

ADPCM didn't get worse. Everyone just stopped trying.  
A small, fast, waveform audio codec that fixes what IMA ADPCM got wrong.  
No transforms. No psychoacoustics. No entropy coding. No excuses.

## what is it?

ADPCM mixed with FLAC in ~450 lines of C.  

## warning

This is work in progress, expect breaking changes.

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
| Marche_Persanne.wav        | 27.31 | 383  | 61.58 dB | 37.60 dB | 1.8026s  | 0.049s   |
| Egyptischer_Marsch.wav     | 32.13 | 623  | 60.17 dB | 37.61 dB | 4.1491s  | 0.090s   |
| suppe_poet_and_peasant.wav | 37.07 | 773  | 58.93 dB | 36.24 dB | 10.3231s | 0.221s   |
| the_four_seasons.wav       | 25.03 | 245  | 62.34 dB | 41.79 dB | 64.2167s | 1.373s   |

**MADPCM** (fast; no lookahead):  
| File                       | RMSE  | MAE  | PSNR     | SNR      | Enc Time |
| -------------------------- | ----- | ---- | -------- | -------- | -------- |
| Marche_Persanne.wav        | 29.02 | 512  | 61.05 dB | 37.07 dB | 0.2169s  |
| Egyptischer_Marsch.wav     | 34.04 | 686  | 59.67 dB | 37.10 dB | 0.5044s  |
| suppe_poet_and_peasant.wav | 39.38 | 800  | 58.40 dB | 35.72 dB | 1.2460s  |
| the_four_seasons.wav       | 26.40 | 247  | 61.88 dB | 41.33 dB | 7.6906s  |

All tested on a Ryzen 5 3600XT running Windows Server 2025 Datacenter build 26100.32230  
Decode speed for MADPCM is the same on slow and fast, this affects only the encoder.  
Encode speed for ADPCM-XQ includes file I/O, tested on a ramdisk to minimalize potential penalty.
File size reduction: ~74% vs WAV.

<details><summary>how a samsung s21 behaves</summary>

```
~/madpcm $ ./test.out encode slow suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
Loaded suppe_poet_and_peasant.wav: 48000 Hz, Stereo
Encoded in 19.806205 seconds (slow mode)
~/madpcm $ ./test.out encode fast suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
Loaded suppe_poet_and_peasant.wav: 48000 Hz, Stereo
Encoded in 2.185003 seconds (fast mode)
~/madpcm $ ./test.out decode suppe_poet_and_peasant_ad.wav suppe_poet_and_peasant_d.wav
Loaded suppe_poet_and_peasant_ad.wav: 2 channels, 48000 Hz, 28272640 samples
Decoded in 0.275335 seconds
```

</details>

## does it sound better?

To me. I think. Better than IMA ADPCM.  
No formal ABX tests were done yet.

## how 2 use it

include `madpcm.c` in your build and include `madpcm.h` in your code.  
preprocessor defines you can use in your project:  
`MADPCM_FREESTANDING` builds the engine for freestanding  
`MADPCM_MEMFUNCS` defines memcpy and memset  
Use with MSVC /Oi- or always define it for freestanding gcc/clang.

## what this does not care about

- standards
- legacy compatibility
- your favorite codec comparison

## roadmap

maybe

## how to contribute

Test it. Report bugs. Open PRs.  
4 spaces, K&R styling.

## q&a

Q: Why make this?  
A: IMA ADPCM works. This works better. IMA's everywhere because it's cheap and simple, not because it's good.

Q: Is this better than Opus/AAC/MP3?  
A: If you squint hard enough and just look at the numbers.

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
A: When it's finished. So soon™ i hope.
