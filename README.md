# Digital Filter Designing - For Noise Reduction

---

## Overview
This project is an in-depth study of bio-signal processing, specifically focusing on filtering and noise reduction techniques applied to ECG signals and brain signals. The report documents the design, implementation, and interpretation of various signal processing methods using MATLAB.

---

## Contents

### 1. Smoothing Filters
One of the most common signal processing tasks is smoothing of the data to reduce high frequency noise arising from electromagnetic interferences, quantization errors and from peripheral physiological signals. Here, the application of moving average filters of order N (MA(N)) and Savitzky-Golay filters of order N and length L’=2L+1 (SG(N,L)) is explored.

- **Moving Average (MA) Filter:**
MA(N) filter can be visualized as a moving window across the signal. A point of the filtered 
signal 𝑦(𝑛) is derived from the average of the windowed points on the input signal 𝑥(𝑛).
This project focuses on,
  - Implementing MA filters and comparing their performance with varying orders.
  - The impact of the filter order on noise reduction and signal distortion(group delay) is analyzed.
![image](https://github.com/user-attachments/assets/0b95e12a-0747-4653-afa8-a809bf84b3d4)

- **Savitzky-Golay (SG) Filter:**
Savitzky-Golay filter fits a polynomial of order 𝑁 to an odd number of data points 𝐿′ = 2𝐿 + 1 
(where 𝐿′ is an odd integer) in a predefined window in a least-squares sense. A unique solution 
requires 𝑁 ≤ 𝐿′ − 1. 
In this project I have covered,
  - Application and parameter optimization(N and L) for better smoothing.
  - Comparative analysis with MA filters.
![image](https://github.com/user-attachments/assets/a99d623a-fbb1-4cfc-ba5d-b585103ad6d6)

### 2. Ensemble Averaging
In the case of overlapping signal and noise spectra, the synchronized averaging technique is an effective and a simple method for noise removal with minimal signal distortions. However, for synchronized averaging to be applicable, there should be input data either having multiple measurements (e.g. EPs) or one signal having repetitive patterns (e.g. ECG). 
This project covers,
- Application of ensemble averaging on noisy ECG signals to improve SNR (This is an example for how to use ensembling average to improve SNR of the signals with repetitive patterns ).
- Use of multiple measurments of Auditory Brain Responses(ABR) to improve SNR.

### 3. FIR Filter Design
Under the context of FIR filter design the project aims,
- Characteristics analysis of different window functions(rectangular, Hanning, Hamming and Blackman ).
- FIR filter implementation using Kaiser windows.
- Design of comb filters to target specific frequencies (50 Hz, 100 Hz, 150 Hz).

### 4. IIR Filters
Similar to FIR filter design project also covers,
- Realization and application of IIR filters.
- Comparison with FIR filters regarding noise suppression, group delay, and frequency response.
- Forward and backward filtering techniques for phase distortion reduction.

---

## Acknowledgments
This project was completed as part of the **BM4152 Bio-Signal Processing** course.
---

