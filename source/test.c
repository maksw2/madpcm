#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <madpcm.h>

typedef struct {
    char riff[4];
    uint32_t fileSize;
    char wave[4];
} RiffHeader;

int load_wav(const char* filename, int16_t** L, int16_t** R, int* rate, int* channels) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        printf("Error: Could not open %s\n", filename);
        return 0;
    }

    RiffHeader h;
    if (fread(&h, sizeof(h), 1, f) != 1 || memcmp(h.riff, "RIFF", 4) != 0 || memcmp(h.wave, "WAVE", 4) != 0) {
        printf("Error: Not a valid WAV file.\n");
        fclose(f);
        return 0;
    }

    char chunkId[4];
    uint32_t chunkSize;
    uint16_t bits = 0, audioFmt = 0, numCh = 0;
    uint32_t dataLen = 0;
    long dataPos = 0;
    int fmtFound = 0;

    while (fread(chunkId, 1, 4, f) == 4) {
        if (fread(&chunkSize, 4, 1, f) != 1) break;
        if (memcmp(chunkId, "fmt ", 4) == 0) {
            fread(&audioFmt, 2, 1, f);
            fread(&numCh, 2, 1, f);
            fread(rate, 4, 1, f);
            fseek(f, 6, SEEK_CUR);
            fread(&bits, 2, 1, f);
            if (chunkSize > 16) fseek(f, chunkSize - 16, SEEK_CUR);
            fmtFound = 1;
        }
        else if (memcmp(chunkId, "data", 4) == 0) {
            dataLen = chunkSize;
            dataPos = ftell(f);
            break;
        }
        else {
            fseek(f, chunkSize, SEEK_CUR);
        }
    }

    if (!fmtFound || dataLen == 0 || bits != 16 || (numCh != 1 && numCh != 2)) {
        printf("Error: Unsupported format (%d ch, %d bits).\n", numCh, bits);
        fclose(f);
        return 0;
    }

    int samples = dataLen / (numCh * 2);
    *channels = numCh;
    *L = malloc(samples * sizeof(int16_t));
    *R = (numCh == 2) ? malloc(samples * sizeof(int16_t)) : NULL;

    fseek(f, dataPos, SEEK_SET);
    for (int i = 0; i < samples; i++) {
        int16_t s;
        fread(&s, 2, 1, f);
        (*L)[i] = s;
        if (numCh == 2) {
            fread(&s, 2, 1, f);
            (*R)[i] = s;
        }
    }

    fclose(f);
    printf("Loaded %s: %d Hz, %s\n", filename, *rate, (numCh == 2) ? "Stereo" : "Mono");
    return samples;
}

int load_wav_madpcm(const char* filename, uint8_t** data, int* size, int* rate, int* totalSamples, int* channels) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        printf("Error: Could not open %s\n", filename);
        return 0;
    }

    RiffHeader h;
    if (fread(&h, sizeof(h), 1, f) != 1 ||
        memcmp(h.riff, "RIFF", 4) != 0 ||
        memcmp(h.wave, "WAVE", 4) != 0) {
        printf("Error: Not a valid RIFF/WAVE file.\n");
        fclose(f);
        return 0;
    }

    char chunkId[4];
    uint32_t chunkSize;
    int fmtFound = 0;
    uint16_t audioFmt = 0;
    uint16_t numChannels = 0;
    uint32_t dataLen = 0;

    while (fread(chunkId, 1, 4, f) == 4) {
        if (fread(&chunkSize, 4, 1, f) != 1) break;

        if (memcmp(chunkId, "fmt ", 4) == 0) {
            fread(&audioFmt, 2, 1, f);
            fread(&numChannels, 2, 1, f);
            fread(rate, 4, 1, f);
            fseek(f, chunkSize - 8, SEEK_CUR);
            fmtFound = 1;
        }
        else if (memcmp(chunkId, "data", 4) == 0) {
            dataLen = chunkSize;
            *data = malloc(dataLen);
            if (fread(*data, 1, dataLen, f) != dataLen) {
                free(*data);
                fclose(f);
                return 0;
            }
            break;
        }
        else {
            fseek(f, chunkSize, SEEK_CUR);
        }
    }

    fclose(f);

    if (!fmtFound || audioFmt != WAVE_FORMAT_MADPCM) {
        printf("Error: File is not in MADPCM format (0x%X).\n", audioFmt);
        return 0;
    }

    *channels = (int)numChannels;
    *size = (int)dataLen;

    // Determine block size based on channel count to calculate total samples
    int blockSize = (numChannels == 2) ? MADPCM_STEREO_SIZE : MADPCM_MONO_SIZE;
    int numBlocks = dataLen / blockSize;
    *totalSamples = numBlocks * MADPCM_BLOCK_SAMPLES;

    printf("Loaded %s: %d channels, %d Hz, %d samples\n", filename, *channels, *rate, *totalSamples);

    return 1;
}

void write_wav_header(FILE* f, int audioFormat, int channels, int rate, int bits, int dataSize) {
    fwrite("RIFF", 1, 4, f);
    uint32_t fileSize = dataSize + 36;
    fwrite(&fileSize, 4, 1, f);
    fwrite("WAVE", 1, 4, f);
    fwrite("fmt ", 1, 4, f);
    uint32_t fmtLen = 16;
    fwrite(&fmtLen, 4, 1, f);
    uint16_t af = audioFormat; fwrite(&af, 2, 1, f);
    uint16_t ch = channels;    fwrite(&ch, 2, 1, f);
    uint32_t sr = rate;        fwrite(&sr, 4, 1, f);
    
    uint16_t blockAlign = (audioFormat == 1) ? (channels * (bits / 8)) : (channels == 2 ? MADPCM_STEREO_SIZE : MADPCM_MONO_SIZE);
    uint32_t byteRate = (audioFormat == 1) ? (rate * blockAlign) : ((rate / MADPCM_BLOCK_SAMPLES) * blockAlign);

    fwrite(&byteRate, 4, 1, f);
    fwrite(&blockAlign, 2, 1, f);
    uint16_t b = bits; fwrite(&b, 2, 1, f);
    fwrite("data", 1, 4, f);
    fwrite(&dataSize, 4, 1, f);
}

void write_wav_pcm(const char* filename, const int16_t* L, const int16_t* R, int samples, int rate) {
    FILE* f = fopen(filename, "wb");
    if (!f) return;

    write_wav_header(f, 1, 2, rate, 16, samples * 4);

    int16_t frame[2];
    for (int i = 0; i < samples; i++) {
        frame[0] = L[i];
        frame[1] = R[i];
        fwrite(frame, sizeof(int16_t), 2, f);
    }
    fclose(f);
    printf("Wrote %s\n", filename);
}

void write_wav_madpcm(const char* filename, const uint8_t* data, int size, int rate) {
    FILE* f = fopen(filename, "wb");
    if (!f) return;

    // Format 0x4D41, 0 bits per sample (compressed)
    write_wav_header(f, WAVE_FORMAT_MADPCM, 2, rate, 0, size);
    fwrite(data, 1, size, f);
    fclose(f);
    printf("Wrote %s (Encoded size: %d bytes)\n", filename, size);
}

// Helper to print usage
void print_usage(const char* exeName) {
    printf("Usage:\n");
    printf("  %s encode [fast|slow] <input.wav> <output.madpcm>\n", exeName);
    printf("  %s decode <input.madpcm> <output.wav>\n", exeName);
}

int main(int argc, char* argv[]) {
    if (argc < 4) { print_usage(argv[0]); return 1; }

    char* mode = argv[1];
    bool isEncode = (strcmp(mode, "encode") == 0);
    bool fast = true;
    char* inputPath, * outputPath;

    if (isEncode) {
        if (argc < 5) { print_usage(argv[0]); return 1; }
        fast = (strcmp(argv[2], "slow") != 0);
        inputPath = argv[3]; outputPath = argv[4];
    }
    else {
        inputPath = argv[2]; outputPath = argv[3];
    }

    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);

    if (isEncode) {
        int16_t* inL, * inR;
        int rate, channels;
        int samples = load_wav(inputPath, &inL, &inR, &rate, &channels);
        if (!samples) return 1;

        int pad = MADPCM_BLOCK_SAMPLES - (samples % MADPCM_BLOCK_SAMPLES);
        if (pad == MADPCM_BLOCK_SAMPLES) pad = 0;
        int totalSamples = samples + pad;

        int16_t* tmpL = realloc(inL, totalSamples * sizeof(int16_t));
        if (tmpL == NULL) {
            free(inL);
            free(inR); // free() on NULL is a no-op
            printf("realloc failed!\n");
            return -1;
        }

        inL = tmpL;
        memset(inL + samples, 0, pad * sizeof(int16_t));
        if (channels == 2) {
            int16_t* tmpR = realloc(inR, totalSamples * sizeof(int16_t));
            if (tmpR == NULL) {
                printf("realloc failed!\n");
                free(inL);
                free(inR);
                return -1;
            }
            inR = tmpR;
            memset(inR + samples, 0, pad * sizeof(int16_t));
        }

        int numBlocks = totalSamples / MADPCM_BLOCK_SAMPLES;
        int blockSize = (channels == 2) ? MADPCM_STEREO_SIZE : MADPCM_MONO_SIZE;
        uint8_t* outBuf = malloc(numBlocks * blockSize);
        if (!outBuf) {
            printf("malloc failed!\n");
            free(inL);
            free(inR);
            return -1;
        }
        int16_t hM[4] = { 0 }, hS[4] = { 0 };

        QueryPerformanceCounter(&start);
        for (int i = 0; i < numBlocks; i++) {
            int off = i * MADPCM_BLOCK_SAMPLES;
            if (channels == 2)
                encodeMADPCMStereoBlock(inL + off, inR + off, outBuf + (i * blockSize), hM, hS, fast);
            else
                encodeMADPCMBlock(inL + off, outBuf + (i * blockSize), hM, fast);
        }
        QueryPerformanceCounter(&end);

        printf("Encoded in %f seconds (%s mode)\n", (double)(end.QuadPart - start.QuadPart) / freq.QuadPart, fast ? "fast" : "slow");
        FILE* f = fopen(outputPath, "wb");
        write_wav_header(f, WAVE_FORMAT_MADPCM, channels, rate, 0, numBlocks * blockSize);
        fwrite(outBuf, 1, numBlocks * blockSize, f);
        fclose(f);

        free(inL); free(inR); free(outBuf);
    }
    else {
        uint8_t* inBuf;
        int inSize, rate, totalSamples, channels;
        if (!load_wav_madpcm(inputPath, &inBuf, &inSize, &rate, &totalSamples, &channels)) return 1;

        int16_t* outL = malloc(totalSamples * sizeof(int16_t));
        if (!outL) {
            printf("malloc failed!\n");
            free(inBuf);
            return -1;
        }
        int16_t* outR = (channels == 2) ? malloc(totalSamples * sizeof(int16_t)) : NULL;
        if (channels == 2 && !outR) {
            printf("malloc failed!\n");
            free(outL);
            free(outR);
            free(inBuf);
            return -1;
        }

        int16_t hM[4] = { 0 }, hS[4] = { 0 };
        int blockSize = (channels == 2) ? MADPCM_STEREO_SIZE : MADPCM_MONO_SIZE;

        QueryPerformanceCounter(&start);
        for (int i = 0; i < totalSamples / MADPCM_BLOCK_SAMPLES; i++) {
            int offS = i * MADPCM_BLOCK_SAMPLES;
            int offB = i * blockSize;
            if (channels == 2) decodeMADPCMStereoBlock(inBuf + offB, outL + offS, outR + offS, hM, hS);
            else decodeMADPCMBlock(inBuf + offB, outL + offS, hM);
        }
        QueryPerformanceCounter(&end);

        printf("Decoded in %f seconds\n", (double)(end.QuadPart - start.QuadPart) / freq.QuadPart);

        FILE* f = fopen(outputPath, "wb");
        write_wav_header(f, 1, channels, rate, 16, totalSamples * channels * 2);
        for (int i = 0; i < totalSamples; i++) {
            fwrite(&outL[i], 2, 1, f);
            if (channels == 2) fwrite(&outR[i], 2, 1, f);
        }
        fclose(f);
        free(outL); if (outR) free(outR); free(inBuf);
    }
    return 0;
}
