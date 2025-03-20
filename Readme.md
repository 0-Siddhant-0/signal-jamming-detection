# DFT Processing for Jammed Signals

Discrete Fourier Transform (DFT) processing techniques to detect signals under jamming conditions. 
The simulation focuses on two challenging scenarios: sum of sinusoids jamming and white Gaussian noise jamming, each requiring different signal processing approaches.

### RF to Discrete-Time Frequency Mapping

The mapping from RF frequency to discrete-time frequency (radians/sample) is defined as:

$$\omega = -\pi + \left(\frac{2\pi}{0.1 \textrm{ GHz}}\right) \times (f_{RF} - 1 \textrm{ GHz})$$

where $\omega$ represents the discrete-time frequency, $f_{RF}$ is the RF frequency with the band ranging from 1.0 GHz ($-\pi$) to 1.1 GHz ($\pi$).

### Complex Signal Generation

The desired signal is generated as a complex exponential:

$$x[n] = e^{j\omega n} = \cos(\omega n) + j\sin(\omega n)$$

where $n$ is the sample index and $\omega$ is the discrete frequency corresponding to one of 18 possible frequencies in the range.

### Jamming Scenarios

#### Scenario I: Sum of Sinusoids

The first jamming scenario uses a sum of 19 evenly spaced sinusoids across the band:

$$x_{jam}[n] = \sum_{k=1}^{19} e^{j\omega_k n}$$

With jamming power set 60 dB higher than the desired signal:

$$x_{jam}[n] = x_{jam}[n] \times \sqrt{10^{60/10}}$$

#### Scenario II: White Gaussian Noise

The white noise jamming is generated with:

$$x_{noise}[n] = \sqrt{\frac{\sigma^2}{2}}(randn[n] + j \cdot randn[n])$$

where $\sigma^2 = 10^{30/10} = 1000$ corresponds to a Signal-to-Noise Ratio (SNR) of -30 dB.

### Signal Detection Techniques

#### Windowing

Three windowing functions were implemented to reduce spectral leakage:
- Hann window: $w[n] = 0.5 \times (1 - \cos(\frac{2\pi n}{N-1}))$
- Hamming window: $w[n] = 0.54 - 0.46 \times \cos(\frac{2\pi n}{N-1})$
- Blackman window: $w[n] = 0.42 - 0.5 \times \cos(\frac{2\pi n}{N-1}) + 0.08 \times \cos(\frac{4\pi n}{N-1})$

#### DFT Averaging

For the white noise scenario, DFT averaging with overlapping windows was applied:

$$X_{avg}[k] = \frac{1}{L} \sum_{i=1}^{L} X_i[k]$$

where $X_i[k]$ represents the DFT of the $i$-th windowed segment with 50% overlap.

## Results

The simulation successfully demonstrates:
- Detection of signals 60 dB below jamming sinusoids using appropriate windowing
- Detection of signals with SNR of -30 dB using DFT averaging
- Comparative analysis of different windowing functions

- `figures/true_and_jam_sig.png`: Visualization of desired signal and jamming signal spectra
- `figures/win_wo.png`: Comparison of detection with and without windowing
- `figures/seg_avg.png`: Demonstration of SNR improvement through segmented DFT averaging
