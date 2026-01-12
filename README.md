# Heat Exchanger Design - Iteration 2

## Overview
This Python script designs a **shell-and-tube heat exchanger** for water vapor condensation. Water vapor entering at 0.0106 MPa and 60°C is fully condensed, transferring 43 kW of heat to cooling water at 33°C inlet temperature.

## Key Features
- **Tube-side calculations**: Reynolds number, Nusselt number, and heat transfer coefficient for water flow
- **Shell-side calculations**: Heat transfer coefficient for vapor condensation using Bonancina correlation
- **Pressure drop analysis**: Calculates friction losses on both tube and shell sides using Kern's method
- **Geometric optimization**: Determines optimal number of passes and tube configuration
- **Visualization**: Generates plots for design parameters and pressure drop relationships

## Main Outputs
- Required heat transfer area: ~2.40 m²
- Tube configuration: 4 passes with optimized tube length
- Overall heat transfer coefficient: ~2050 W/(m²·K)
- Pressure drop predictions for both sides
- Summary table of all calculated parameters

## Dependencies
- `numpy>=1.20.0`
- `pandas>=1.3.0`
- `matplotlib>=3.4.0`
- `scipy>=1.7.0`

## Installation
Install required packages using:
```bash
pip install -r requirements.txt
```

## Usage
Run the script directly:
```bash
python Iteration_2.py
```

The script will print a comprehensive summary table and generate visualization plots for design validation.

## Design Approach
1. Calculates LMTD (Log Mean Temperature Difference)
2. Determines required heat transfer area
3. Sizes tubes based on water velocity and pressure drop constraints
4. Iterates on shell diameter and baffle spacing
5. Validates heat transfer coefficients and pressure drops
6. Provides visualization of key design relationships

## Author
Kaleb Kassa
