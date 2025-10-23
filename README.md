# Look-Ahead IIR at 125 MHz on NI FlexRIO (Kintex-7)

A k=1 look-ahead biquad that sustains a 125 MHz SCTL on PXIe-7972R + NI-5782R, preserves dynamics (only a pure +1-sample delay), and achieves ≈ **−81.8 dB NMSE** vs a fixed-point MATLAB reference. Timeseries stream via DMA; host plots/validates.

**Hardware/Tools**
- NI PXIe-7972R (Kintex-7) + NI-5782R I/O, I/O Module Clock at 125 MHz (SCTL).
- LabVIEW FPGA (2023+), NI drivers; MATLAB (for verification).

**Design highlights**
- All multiplies map to DSP48E1; balanced adder tree; uniform 1-tick pipeline on every branch ⇒ effect is a pure *z⁻¹* delay. Signals carried as Q2.13 end-to-end; single final cast.
- High-rate streams via Target→Host DMA; host builds waveforms with dt = 8 ns × Ndec.

**Quickstart**
1) Open `labview/LookAheadIIR.lvproj`, bind hardware, build or bind the `.lvbitx`.
2) Run `host` VI → you should see Raw/AI/Filtered streams.
3) `matlab/verify_iir.m` to reproduce parity, PSD overlap, and error PSD with the fixed-point reference (no gain/offset fit; accounts for the known +1 sample).

**Key results**
- Relative RMS error ≈ 8.13×10⁻⁵, NMSE ≈ −81.8 dB, PSD overlap 1–60 MHz; added dynamics: only +1 sample.

**Reference**
See the full reproducible report (UNSW practicum) in `/docs` for the exact DMA names, numerics, and timing-closure checklist.


