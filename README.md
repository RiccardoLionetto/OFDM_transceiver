# OFDM Transceiver

This repository contains a Matlab implementation designed to communicate using Orthogonal Frequency-Division Multiplexing (OFDM).

### Introduction to OFDM
OFDM is a modulation technique widely used in digital communication systems. It allows to transmit data over dispersive channels while mitigating the effects of interference and multipath propagation.

### Implemented Features
Main features are:
- continuous phase tracking (Viterbi-Viterbi algorithm)
- channel estimation and correction
- data pseudo-randomization
- data windowing
- pilot pattern: block type

The above techniques allowed to obtain robustness against PAPR problem, carrier frequency offset, multipath phenomenon and channel fading.

### Hardware Compatibility
The algorithm is designed to be compatible with simple and readily available hardware components, requiring only a microphone and a pair of speakers. As such, acoustic waves were utilized to test the code.

On one side this choice bring accessibility and cost-effectiveness, on the other it imposes to be careful on the limited frequency range of speakers and microphone. Moreover, due to the short wavelenghts involved, little movements of the hardware during transmission can introduce phase shifts and distortions.

### Results
- <b>Viterbi-Viterbi phase correction</b>: to highlight the positive contribution of this feature a carrier frequency offset was artificially forced in the received signal.
![Viterbi-Viterbi correction](/img/ViterbiViterbi_Correction.png)

