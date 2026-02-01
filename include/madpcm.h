#pragma once

#include <stdint.h>

// Core MADPCM functions
int encodeMADPCMBlock(const int16_t* inSamples, uint8_t* outBuffer, int16_t* history);
void decodeMADPCMBlock(const uint8_t* inBuffer, int16_t* outSamples, int16_t* history);
int encodeMADPCMStereoBlock(const int16_t* inL, const int16_t* inR, uint8_t* outBuffer, int16_t* historyM, int16_t* historyS);
void decodeMADPCMStereoBlock(const uint8_t* inBuffer, int16_t* outL, int16_t* outR, int16_t* historyM, int16_t* historyS);

#define MADPCM_BLOCK_SAMPLES 1024
#define MADPCM_HEADER_SIZE 16 
#define MADPCM_PAYLOAD_SIZE 512
// Mono Block = 528 bytes
#define MADPCM_MONO_SIZE (MADPCM_HEADER_SIZE + MADPCM_PAYLOAD_SIZE) 
// Stereo Block = 1056 bytes
#define MADPCM_STEREO_SIZE (MADPCM_MONO_SIZE * 2)
#define WAVE_FORMAT_MADPCM 0x4D41 // 'MA' in hex
