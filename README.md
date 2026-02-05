# Uranium-238 Decay Chain Analysis

![Language](https://img.shields.io/badge/language-Python%20%7C%20Fortran-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-Completed-success)

## ☢️ Project Overview

This project provides a comprehensive computational analysis of the radioactive decay chain starting from **Uranium-238** and ending at the stable isotope **Lead-206**. The study simulates the isotopic evolution over the age of the Earth (approx. 4.503 billion years).

The primary goal is to compare **Analytical** methods (exact solutions) with **Numerical** integration methods to demonstrate numerical stability challenges in "stiff" differential equations caused by vastly different half-lives.

### 🔬 Key Features
* **Analytical Solution:** Implementation of **Bateman Equations** for exact mass calculations.
* **Numerical Analysis:**
  * **Forward Euler (Explicit):** Demonstrates stability limitations when time steps ($\Delta t$) are large.
  * **Implicit Euler (Fortran):** Provides a robust solution for stiff decay chains.
* **Reverse Geochronology:** Calculates the initial U-238 mass at Earth's formation based on present-day Pb-206 abundance.
* **High Precision:** Uses `numpy.longdouble` in Python and `double precision` in Fortran to minimize floating-point errors.

---

## 📊 Visuals & Results

The analysis includes logarithmic plots tracking the mass evolution of U-238, U-234, Th-230, Ra-226, Pb-210, and Pb-206.

![Decay Chain Graph](docs/graph_preview.png)

The study highlights that while Uranium-238 decays slowly, intermediate isotopes like **Pb-210** (half-life approx. 22 years) require extremely small time steps or implicit methods for accurate numerical simulation.

---

## 📂 Project Structure

    uranium-decay-chain-analysis/
    ├── docs/
    │   ├── U238DecayChain.pdf       # Full Project Report
    │   └── graph_preview.png        # Preview image for README
    ├── src/
    │   ├── python/
    │   │   └── decay_chain_analytic.py  # Bateman & Forward Euler implementation
    │   └── fortran/
    │       └── implicit_euler.f90       # Implicit Euler implementation
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

The decay chain is modeled using a system of first-order differential equations.


The general equation for the decay rate is:

$$\frac{dN_i}{dt} = \lambda_{i-1}N_{i-1} - \lambda_i N_i$$

Where:
* $\lambda$ is the decay constant ($\ln(2) / t_{1/2}$).
* $N$ is the number of atoms.

**The Decay Chain:**
`U-238` -> `U-234` -> `Th-230` -> `Ra-226` -> `Pb-210` -> `Pb-206 (Stable)`

---

## 👨‍💻 Author

**Emre Sakarya**
* Hacettepe University, Department of Nuclear Engineering
* Project: NEM 393 Engineering Project II

---

*For detailed physics, derivations, and error analysis, please refer to the [Project Report](docs/U238DecayChain.pdf).*
