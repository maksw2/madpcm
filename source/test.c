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

int load_wav_stereo(const char* filename, int16_t** L, int16_t** R, int* rate) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        printf("Error: Could not open %s\n", filename);
        return 0;
    }

    RiffHeader h;
    if (fread(&h, sizeof(h), 1, f) != 1 ||
        memcmp(h.riff, "RIFF", 4) != 0 ||
        memcmp(h.wave, "WAVE", 4) != 0) {
        printf("Error: Not a valid WAV file.\n");
        fclose(f);
        return 0;
    }

    char chunkId[4];
    uint32_t chunkSize;
    uint16_t channels = 0, bits = 0, audioFmt = 0;
    uint32_t dataLen = 0;
    long dataPos = 0;
    int fmtFound = 0;

    while (fread(chunkId, 1, 4, f) == 4) {
        if (fread(&chunkSize, 4, 1, f) != 1) break;

        if (memcmp(chunkId, "fmt ", 4) == 0) {
            fread(&audioFmt, 2, 1, f);
            fread(&channels, 2, 1, f);
            fread(rate, 4, 1, f);
            fseek(f, 6, SEEK_CUR);
            fread(&bits, 2, 1, f);

            int readSoFar = 16;
            if (chunkSize > (uint32_t)readSoFar) fseek(f, chunkSize - readSoFar, SEEK_CUR);
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

    if (!fmtFound || dataLen == 0) {
        printf("Error: Could not find 'fmt ' or 'data' chunks.\n");
        fclose(f);
        return 0;
    }

    if (channels != 2 || bits != 16) {
        printf("Error: Input must be 16-bit Stereo PCM. Found: %d ch, %d bits.\n", channels, bits);
        fclose(f);
        return 0;
    }

    int samples = dataLen / 4;
    *L = malloc(samples * sizeof(int16_t));
    *R = malloc(samples * sizeof(int16_t));

    fseek(f, dataPos, SEEK_SET);
    int16_t frame[2];
    for (int i = 0; i < samples; i++) {
        if (fread(frame, sizeof(int16_t), 2, f) < 2) break;
        (*L)[i] = frame[0];
        (*R)[i] = frame[1];
    }

    fclose(f);

    // Normalized Output: 1. Channels, 2. Size (bytes), 3. Samples
    printf("Loaded %s: Channels: %d, Size: %u bytes, Samples: %d (%d Hz)\n",
        filename, channels, dataLen, samples, *rate);

    return samples;
}

int load_wav_madpcm(const char* filename, uint8_t** data, int* size, int* rate, int* totalSamples) {
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
    uint32_t dataLen = 0;

    while (fread(chunkId, 1, 4, f) == 4) {
        if (fread(&chunkSize, 4, 1, f) != 1) break;

        if (memcmp(chunkId, "fmt ", 4) == 0) {
            fread(&audioFmt, 2, 1, f);
            fseek(f, 2, SEEK_CUR); // Skip channels
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

    *size = (int)dataLen;
    int numBlocks = dataLen / MADPCM_STEREO_SIZE;
    *totalSamples = numBlocks * MADPCM_BLOCK_SAMPLES;

    // Normalized Output: 1. Channels, 2. Size (bytes), 3. Samples
    printf("Loaded %s: Channels: 2, Size: %d bytes, Samples: %d (%d Hz)\n",
        filename, *size, *totalSamples, *rate);

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
    
    // Block Align / Byte Rate
    uint32_t byteRate;
    uint16_t blockAlign;

    if (audioFormat == 1) { // PCM
        blockAlign = channels * (bits / 8);
        byteRate = rate * blockAlign;
    } else { // MADPCM (0x4D41)
        // 1024 samples -> 1056 bytes stereo
        // Effective byte rate
        blockAlign = MADPCM_STEREO_SIZE;
        byteRate = (rate / MADPCM_BLOCK_SAMPLES) * MADPCM_STEREO_SIZE;
    }

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
    printf("  %s encode <input.wav> <output.madpcm>\n", exeName);
    printf("  %s decode <input.madpcm> <output.wav>\n", exeName);
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    char* mode = argv[1];
    char* inputPath = argv[2];
    char* outputPath = argv[3];

    LARGE_INTEGER frequency, start, end;
    QueryPerformanceFrequency(&frequency);

    if (strcmp(mode, "encode") == 0) {
        int16_t *inL, *inR;
        int rate;
        
        int samples = load_wav_stereo(inputPath, &inL, &inR, &rate);
        if (!samples) {
            fprintf(stderr, "Failed to load input WAV: %s\n", inputPath);
            return 1;
        }

        // Pad to multiple of MADPCM_BLOCK_SAMPLES
        int pad = MADPCM_BLOCK_SAMPLES - (samples % MADPCM_BLOCK_SAMPLES);
        if (pad == MADPCM_BLOCK_SAMPLES) pad = 0;
        int totalSamples = samples + pad;

        if (pad > 0) {
            int16_t* tmpL = realloc(inL, totalSamples * sizeof(int16_t));
            int16_t* tmpR = realloc(inR, totalSamples * sizeof(int16_t));

            if (tmpL == NULL || tmpR == NULL) {
                printf("realloc failed!\n");
                // free() on null is a guaranteed no-op
                free(tmpL);
                free(tmpR);
                return -1;
            }

            inL = tmpL;
            inR = tmpR;

            memset(inL + samples, 0, pad * sizeof(int16_t));
            memset(inR + samples, 0, pad * sizeof(int16_t));
        }

        int numBlocks = totalSamples / MADPCM_BLOCK_SAMPLES;
        // Stereo Block = (Header(10) + Payload(512)) * 2 = 1056 bytes
        int maxEncodedSize = numBlocks * MADPCM_STEREO_SIZE;
        uint8_t* encodedBuf = malloc(maxEncodedSize);
        int encodedSize = 0;
        int16_t histM[4] = {0}, histS[4] = {0};

        printf("Encoding %d blocks...\n", numBlocks);
        QueryPerformanceCounter(&start);

        for (int i = 0; i < numBlocks; i++) {
            int offset = i * MADPCM_BLOCK_SAMPLES;
            int bytes = encodeMADPCMStereoBlock(
                inL + offset, 
                inR + offset, 
                encodedBuf + encodedSize, 
                histM, histS
            );
            encodedSize += bytes;
        }

        QueryPerformanceCounter(&end);
        printf("Encoding took %f seconds\n", (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart);

        write_wav_madpcm(outputPath, encodedBuf, encodedSize, rate);
        
        free(inL); free(inR); free(encodedBuf);
    } 
    else if (strcmp(mode, "decode") == 0) {
        uint8_t* encodedBuf;
        int encodedSize, rate, totalSamples;

        if (!load_wav_madpcm(inputPath, &encodedBuf, &encodedSize, &rate, &totalSamples)) {
            fprintf(stderr, "Failed to load MADPCM file: %s\n", inputPath);
            return 1;
        }

        int numBlocks = totalSamples / MADPCM_BLOCK_SAMPLES;
        int16_t* outL = malloc(totalSamples * sizeof(int16_t));
        int16_t* outR = malloc(totalSamples * sizeof(int16_t));
        int16_t histM[4] = {0}, histS[4] = {0};

        printf("Decoding %d blocks...\n", numBlocks);
        QueryPerformanceCounter(&start);

        int readPtr = 0;
        for (int i = 0; i < numBlocks; i++) {
            decodeMADPCMStereoBlock(
                encodedBuf + readPtr, 
                outL + (i * MADPCM_BLOCK_SAMPLES), 
                outR + (i * MADPCM_BLOCK_SAMPLES), 
                histM, histS
            );
            readPtr += MADPCM_STEREO_SIZE; // Fixed size for stereo blocks
        }

        QueryPerformanceCounter(&end);
        printf("Decoding took %f seconds\n", (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart);

        write_wav_pcm(outputPath, outL, outR, totalSamples, rate);

        free(outL); free(outR); free(encodedBuf);
    } 
    else {
        print_usage(argv[0]);
        return 1;
    }

    return 0;
}
