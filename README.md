# 🌌 Astrophysical Simulations by @me

This repository contains two astrophysical simulation projects developed as part of a research study in computational astrophysics:

1. **Iron Diffusion in the Intracluster Medium (ICM)** – Perseus Cluster modeling
2. **Supernova Remnant (SNR) Evolution** – Explosion and shock propagation in the ISM

---

## 📘 Project 1: Iron Diffusion in the ICM

### 🧠 Overview

Simulates the diffusion of iron in the ICM of the Perseus galaxy cluster by solving hydrostatic equilibrium, dark matter potential, and diffusion equations.

### 🔧 Features

- 🌀 NFW dark matter profile
- 🧪 Hydrostatic equilibrium for gas density profile
- 🌡️ Variable temperature from X-ray observations
- 🧲 Iron diffusion with/without SN source terms
- 📊 Mass conservation and diffusion diagnostics
- 🔁 Bisection method for baryonic fraction tuning

### 🧪 Simulation Modes

- **Start with observed Fe abundance** – studies diffusion alone
- **Start with zero Fe and add SN injection** – models iron enrichment

### 🖥️ User Inputs

- Total evolution time (Gyr)
- Output interval (Gyr)
- Initial iron abundance option
- Supernova rate (SN per century)

### 📈 Output

- Iron abundance over time (`Z_Fe(r, t)`)
- Mass of Fe in different radial shells
- Diffusion timescale estimations

---

## 🌠 Project 2: Supernova Remnant (SNR) Evolution

### 🧠 Overview

Simulates the evolution of a supernova remnant (SNR) in the ISM using hydrodynamic equations. Compares results against the analytical Sedov-Taylor solution.

### 🔧 Features

- 📏 1D radial grid (500–5000 points)
- 💥 Thermal energy injection (E = 10^51 erg)
- 🔬 Energy partitioning: kinetic vs thermal
- 🌌 Radiative cooling (optional)
- 🌟 X-ray luminosity calculation

### 🧪 Simulation Modes

- **Adiabatic expansion** – classic Sedov phase
- **With radiative cooling** – realistic late-stage SNR
- **Varying ISM conditions** – tests physical sensitivity

### 🖥️ User Inputs

- Grid resolution
- ISM temperature and density
- Cooling switch
- Simulation time and outputs

### 📈 Output

- Shock radius & velocity vs time
- Density/temperature/velocity profiles
- Energy conservation validation
- X-ray luminosity evolution

---

## ⚙️ Requirements

- Fortran compiler (e.g., `gfortran`)
- Python or gnuplot for post-processing plots
- Memory considerations for high-res 