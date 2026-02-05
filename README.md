# Uranium-238 Decay Chain Analysis

![Language](https://img.shields.io/badge/language-Python%20%7C%20Fortran-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-Completed-success)

## ☢️ Project Overview

[cite_start]This project provides a comprehensive computational analysis of the radioactive decay chain starting from **Uranium-238** and ending at the stable isotope **Lead-206**[cite: 25]. [cite_start]The study simulates the isotopic evolution over the age of the Earth (approx. 4.503 billion years)[cite: 40].

[cite_start]The primary goal is to compare **Analytical** methods (exact solutions) with **Numerical** integration methods to demonstrate numerical stability challenges in "stiff" differential equations caused by vastly different half-lives[cite: 32].

### 🔬 Key Features
* [cite_start]**Analytical Solution:** Implementation of **Bateman Equations** for exact mass calculations[cite: 30].
* **Numerical Analysis:**
  * [cite_start]**Forward Euler (Explicit):** Demonstrates stability limitations when time steps ($\Delta t$) are large[cite: 17].
  * [cite_start]**Implicit Euler (Fortran):** Provides a robust solution for stiff decay chains[cite: 17].
* [cite_start]**Reverse Geochronology:** Calculates the initial U-238 mass at Earth's formation based on present-day Pb-206 abundance[cite: 15].
* **High Precision:** Uses `numpy.longdouble` in Python and `double precision` in Fortran to minimize floating-point errors.

---

## 📊 Visuals & Results

[cite_start]The analysis includes logarithmic plots tracking the mass evolution of U-238, U-234, Th-230, Ra-226, Pb-210, and Pb-206[cite: 74].

[cite_start]The study highlights that while Uranium-238 decays slowly, intermediate isotopes like **Pb-210** (half-life approx. 22 years) require extremely small time steps or implicit methods for accurate numerical simulation[cite: 26, 237].

---

## 📂 Project Structure

    uranium-decay-chain-analysis/
    ├─── U238DecayChain.pdf       # Full Project Report & detailed graphs
    │ 
    ├── src/
    │   ├── python/
    │   │   └── decay_chain_analytic.py  # Bateman & Forward Euler implementation
    │   └── fortran/
    │       └── implicit_euler.f90       # Implicit Euler implementation (High Stability)
    └── README.md

---

## 🚀 How to Run

### Prerequisites
* **Python:** NumPy is required.
* **Fortran:** GFortran or any standard Fortran compiler.

### 1. Running the Python Simulation
This script runs the Bateman analytical solution and the Forward Euler numerical test.

    pip install numpy
    python src/python/decay_chain_analytic.py

### 2. Running the Fortran Simulation
This script runs the Implicit Euler method, which is stable even at large time steps (e.g., dt = 1,000,000 years).

    cd src/fortran
    gfortran implicit_euler.f90 -o decay_sim
    ./decay_sim

---

## 🧮 Mathematical Background

[cite_start]The decay chain is modeled using a system of first-order differential equations[cite: 52]:

$$\frac{dN_i}{dt} = \lambda_{i-1}N_{i-1} - \lambda_i N_i$$

Where:
* [cite_start]$\lambda$ is the decay constant ($\ln(2) / t_{1/2}$)[cite: 48].
* $N$ is the number of atoms.

**The Decay Chain:**
[cite_start]`U-238` -> `U-234` -> `Th-230` -> `Ra-226` -> `Pb-210` -> `Pb-206 (Stable)`[cite: 13].

---

## 👨‍💻 Author

[cite_start]**Emre Sakarya** [cite: 8]
* [cite_start]Hacettepe University, Department of Nuclear Engineering [cite: 1, 2]
* [cite_start]Project: NEM 393 Engineering Project II [cite: 3]

---

*For detailed physics, derivations, and error analysis, please refer to the [Project Report](docs/U238DecayChain.pdf).*
