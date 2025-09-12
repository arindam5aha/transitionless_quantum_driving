## Transitionless Phonon Assisted Quantum Driving

This repository is related to the study of transitionless quantum driving, including phonon-assisted processes. For a detailed theoretical background and results, see the included PDF:

- **Transitionless Phonon Assisted Quantum Driving** (`Transitionless Phonon assited.pdf`)

The document provides:
- An overview of transitionless quantum driving techniques
- The role of phonons in quantum state transfer
- Mathematical models and results
- References to experimental and theoretical literature

**How to use:**
Read the PDF for a comprehensive understanding of the physical principles and methods implemented in the code files of this repository. The PDF serves as a scientific foundation for the simulations and analyses performed here.

## Code Files

The following Python scripts are included for simulation and analysis:

- `Fidelity TQD.py`: Calculates fidelity for transitionless quantum driving protocols.
- `g1 g2.py`: Defines and analyzes coupling parameters.
- `Population Adiabatic.py`: Simulates population dynamics under adiabatic conditions.
- `Population TQD.py`: Simulates population dynamics using transitionless quantum driving.
- `TQD Fidelity Map.py`: Generates fidelity maps for different parameter regimes.

## QuTiP Reference

This project uses [QuTiP](http://qutip.org/), the Quantum Toolbox in Python, for simulating quantum systems. QuTiP provides efficient numerical routines for solving the dynamics of open and closed quantum systems.

**To install QuTiP:**
```bash
pip install qutip
```

Refer to the QuTiP documentation for usage details: https://qutip.org/docs/latest/