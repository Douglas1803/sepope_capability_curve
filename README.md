Autor: Douglas Ficagna
E-mail: douglas.elf@outlook.com
LinkedIn: https://www.linkedin.com/in/douglas-ficagna/

# Synchronous Generator Capability Curve Calculation  
## Generic Hydroelectric Generator Model with Excitation Limiters

This repository presents a computational methodology in Python for constructing the **capability curve (P-Q diagram)** of a synchronous generator considering not only classical machine limits, but also **realistic excitation system constraints**.

The project was inspired by academic studies involving hydroelectric power plant operation and aims to provide a **generic, reproducible, and public-safe model** for educational, research, and engineering applications.

---

# Objective

Develop a computational tool capable of representing the operational capability region of a synchronous generator by incorporating:

- Stator current limit (SCL)
- Overexcitation limiter (OEL)
- Minimum excitation limiter (MEL)
- Underexcitation limiter (UEL)
- Load angle limiter
- Magnetic saturation effects
- Practical operating constraints in the P-Q plane

---

# Motivation

Traditional capability curves often represent only theoretical thermal and stability boundaries.  
In real-world operation, especially under reduced dispatch and high reactive support requirements, excitation system limiters frequently define the true operational boundary.

This project provides a more realistic representation of synchronous generator capability by combining:

**Electrical Parameters + Saturation + Excitation Limiters → Practical Capability Curve**

---

# Generator Model

This repository uses **generic synchronous generator parameters** based on typical literature ranges (e.g., Kundur, synchronous machine modeling references), avoiding proprietary or confidential operational data.

### Base Values:
- Rated Power: **100 MVA**
- Terminal Voltage: **13.8 kV**
- Frequency: **60 Hz**
- Power Factor: **0.95**

### Example Parameters:
- Xd = 1.05 pu
- Xq = 0.70 pu
- Xd’ = 0.30 pu
- Xd’’ = 0.20 pu
- Ra = 0.003 pu

---

# Main Features

## Capability Curve Components:
- Classical stator thermal limit
- Rotor field current limit
- Saturation-adjusted field model
- OEL curve
- MEL curve
- UEL representation
- Load angle restriction
- Operational validation points

---

# Repository Structure

```bash
project/
│
├── dados.py              # Generator and excitation system data
├── main.py               # Main execution script
├── curvas.py             # Calculation of curves
├── plotagem.py           # Curve plot
├── calc_iterativo.py     # Calculation of field current limits
├── calc_ANG_CARGA.py     # Calculation of load angle and subexcitation limiters
├── requirements.txt      # Required Python libraries
└── README.md             # Project documentation
