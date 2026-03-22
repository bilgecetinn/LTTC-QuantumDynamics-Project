# LTTC-QuantumDynamics-Project
Summary: Fortran 90 implementation of quantum wave packet propagation using the split-operator method. Built as a graded assignment for Prof. Arjan Berger's Quantum Dynamics module at LTTC Winter School January 2025 (TCCM Master's). Extends tutorial code with harmonic oscillator eigenstates, Hermite polynomials, and superposition analysis.

Quantum Dynamics: Harmonic Oscillator Wave Packet Propagation
LTTC Winter School — January 2025 | TCCM Master's Programme
Aspet, Toulouse, France
Assignment submitted to: Prof. Arjan Berger (LCPQ, Université Paul Sabatier)

Overview
This project extends a quantum dynamics tutorial developed by Prof. Arjan Berger for the LTTC Winter School. The original framework propagates a Gaussian wave packet using the split-operator method. My contribution replaces the Gaussian initial state with a wave packet constructed as a superposition of quantum harmonic oscillator eigenstates, and investigates how different coefficient choices affect the dynamics.
The full assignment report is included: LTTC-ArjanBerger-QDProject-BilgeCetin-HW-report.pdf

Background & Attribution
The base propagation code and tutorial were created by Prof. Arjan Berger, inspired by:

Tanner, J. J., J. Chem. Education 67, 917 (1990)
Original Fortran program by Anders Sandvik (Boston University)

The tutorial exercises covered Gaussian wave packet propagation, norm and energy conservation, autocorrelation functions, and eigenvalue spectra. This repository contains my graded homework extension built on top of that framework.

My Contributions
Task 1 — Hermite Polynomial Recursion
Implemented a hermite_poly subroutine computing Hermite polynomials Hₙ(y) via the three-term recursion:
H₀(y) = 1
H₁(y) = 2y
Hₙ(y) = 2y·Hₙ₋₁(y) − 2(n−1)·Hₙ₋₂(y)
Task 2 — Harmonic Oscillator Eigenstates
Modified initpsi to compute eigenstates φₙ(x) for n = 0 to 4 using the analytic formula:
φₙ(x) = 1/√(2ⁿ n!) · (mω/πℏ)^(1/4) · exp(−mωx²/2ℏ) · Hₙ(√(mω/ℏ) · x)
Task 3 — Superposition Wave Packet
Built the initial wavefunction as a weighted sum:
Ψ₀(x) = Σ c(n) · φₙ(x),  n = 0..4
Coefficients are read from the external coefficients file, making it easy to test different configurations without recompiling.
Task 4 — Normalization
Added a normalization block at the end of initpsi that computes the norm via discrete integration and rescales Ψ₀(x) to satisfy ∫|Ψ₀(x)|²dx = 1.
Tasks 5 & 6 — Physical Validation & Case Studies
Validated the code by propagating a pure eigenstate (only one non-zero coefficient) and confirming that the probability density remains unchanged over time — only a phase factor is acquired, as expected.
Three coefficient configurations were then compared:
CaseCoefficientsSymmetryBehaviour1All equal to 1NoneComplex beating pattern, multiple interference peaks2Even only (c₀, c₂, c₄ = 1)Even (symmetric about x=0)Regular, symmetric oscillation3Odd only (c₁, c₃ = 1)Odd (antisymmetric)Node at x=0, period set by E₃−E₁ gap

Repository Contents
FileDescriptionpropagate-modified.f90Main propagation program (my modified version)analysis.f90Extended analysis: norm, energy, autocorrelation, eigenvalue spectrumgraphics.f90Snapshot generation subroutinesdfft.fFast Fourier Transform auxiliary routinescoefficientsInput file: superposition coefficients c(n)potentialInput file: potential type and parametersLTTC-ArjanBerger-QDProject-BilgeCetin-HW-report.pdfFull submitted assignment report

How to Compile and Run
Requires gfortran and optionally ImageMagick for animated GIF output.
bash# Compile
gfortran propagate-modified.f90 graphics.f90 dfft.f -o propagate

# Run (from a working directory containing the input files)
./propagate

# Convert snapshots to animated GIF (optional)
convert -set dispose previous -delay 20 psi*.ps output.gif

Key Concepts

Quantum wave packet propagation
Split-operator method (Fourier-space kinetic operator)
Hermite polynomials and harmonic oscillator eigenstates
Superposition states and interference
Wavefunction normalization
Eigenvalue spectra via autocorrelation function


Academic Context
This work was completed as part of the Quantum Dynamics module at the LTTC Winter School, January 2025, held in Aspet, France — an intensive course within the TCCM European Master's programme (Theoretical Chemistry and Computational Modelling).
Submitted by: Bilge Emek Cetin
