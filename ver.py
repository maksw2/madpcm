import numpy as np
import wave
import sys
import os

def calculate_metrics_optimized(file_orig, file_proc):
    # Chunk size: 100MB of raw 16-bit data per file
    # (Results in ~800MB of float64 data during processing)
    CHUNK_BYTES = 100 * 1024 * 1024 

    with wave.open(file_orig, 'rb') as w_orig, wave.open(file_proc, 'rb') as w_proc:
        params = w_orig.getparams()
        n_frames = min(w_orig.getnframes(), w_proc.getnframes())
        sample_width = w_orig.getsampwidth()
        
        # Calculate frames per chunk based on 2 bytes per sample (int16)
        # and number of channels
        bytes_per_frame = sample_width * params.nchannels
        frames_per_chunk = CHUNK_BYTES // bytes_per_frame

        sum_sq_diff = 0.0
        sum_sq_signal = 0.0
        max_abs_diff = 0.0
        total_samples = 0

        print(f"Processing {n_frames / params.framerate:.2f} seconds of audio...")

        frames_read = 0
        while frames_read < n_frames:
            to_read = min(frames_per_chunk, n_frames - frames_read)
            
            # 1. Read raw bytes
            buf_o = w_orig.readframes(to_read)
            buf_p = w_proc.readframes(to_read)
            
            # 2. Convert to float64 for high-precision math
            # This is where the memory expands 4x
            orig_chunk = np.frombuffer(buf_o, dtype=np.int16).astype(np.float64)
            proc_chunk = np.frombuffer(buf_p, dtype=np.int16).astype(np.float64)

            # 3. Vectorized Math (Numpy uses SIMD here for speed)
            diff = orig_chunk - proc_chunk
            
            sum_sq_diff += np.sum(diff ** 2)
            sum_sq_signal += np.sum(orig_chunk ** 2)
            
            current_max = np.max(np.abs(diff))
            if current_max > max_abs_diff:
                max_abs_diff = current_max
            
            total_samples += len(orig_chunk)
            frames_read += to_read
            
            # Progress indicator
            progress = (frames_read / n_frames) * 100
            print(f"\rProgress: {progress:.1f}%", end="", flush=True)

        print("\nFinalizing calculations...")

        # Final Metrics
        mse = sum_sq_diff / total_samples
        rmse = np.sqrt(mse)
        psnr = 20 * np.log10(32768.0 / rmse) if rmse > 0 else float('inf')
        snr = 10 * np.log10(sum_sq_signal / sum_sq_diff) if sum_sq_diff > 0 else float('inf')

        return psnr, rmse, snr, max_abs_diff

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python ver.py original.wav decompressed.wav")
        sys.exit(1)

    try:
        psnr, rmse, snr, mae = calculate_metrics_optimized(sys.argv[1], sys.argv[2])

        print(f"RMSE (Avg Error):  {rmse:.4f}")
        print(f"MAE (Peak Error):  {mae:.0f}")
        print(f"PSNR:              {psnr:.2f} dB")
        print(f"SNR:               {snr:.2f} dB")
    except Exception as e:
        print(f"\nError: {e}")