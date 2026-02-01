import numpy as np
import wave
import sys

def load_wav(filepath):
    with wave.open(filepath, 'rb') as w:
        params = w.getparams()
        frames = w.readframes(params.nframes)
        # Convert buffer to int16 and cast to float64 for precision
        samples = np.frombuffer(frames, dtype=np.int16).astype(np.float64)
        return samples, params.framerate

def calculate_metrics(original, processed):
    # Ensure they are the same length
    min_len = min(len(original), len(processed))
    original = original[:min_len]
    processed = processed[:min_len]

    # Calculate differences
    diff = original - processed
    
    # 1. MAE (Maximum Absolute Error)
    mae = np.max(np.abs(diff))

    # 2. MSE & RMSE
    mse = np.mean(diff ** 2)
    if mse == 0:
        return float('inf'), 0, float('inf'), 0
    
    rmse = np.sqrt(mse)
    
    # 3. PSNR (Peak Signal-to-Noise Ratio)
    # Using 32768 as the peak for 16-bit signed audio
    max_val = 32768.0
    psnr = 20 * np.log10(max_val / rmse)

    # 4. SNR (Signal-to-Noise Ratio)
    # Ratio of signal power to noise power
    signal_power = np.mean(original ** 2)
    if signal_power == 0:
        snr = 0 # Avoid division by zero for silent files
    else:
        snr = 10 * np.log10(signal_power / mse)
    
    return psnr, rmse, snr, mae

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python ver.py original.wav decompressed.wav")
        sys.exit(1)

    orig_s, rate = load_wav(sys.argv[1])
    proc_s, _ = load_wav(sys.argv[2])

    psnr, rmse, snr, mae = calculate_metrics(orig_s, proc_s)

    print(f"--- Codec Verification Report ---")
    print(f"RMSE (Average Error): {rmse:.4f} units")
    print(f"MAE (Max Peak Error): {mae:.0f} units")
    print(f"PSNR (Fixed Peak):    {psnr:.2f} dB")
    print(f"SNR (Signal-Based):   {snr:.2f} dB")
    