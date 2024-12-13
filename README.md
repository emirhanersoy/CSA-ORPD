# Crystal Structure Algorithm for IEEE 14-Bus ORPD

This project addresses the Optimal Reactive Power Dispatch (ORPD) problem in the IEEE 14-Bus system using the Crystal Structure Algorithm (CSA). The goal is to minimize power losses while adhering to system constraints.

## Description

The Crystal Structure Algorithm is implemented to optimize reactive power dispatch in power systems. This implementation specifically targets the IEEE 14-Bus test system, providing a solution for minimizing system losses while maintaining operational constraints.

## Project Structure
├── CSA.py                     # Main algorithm implementation

├── IEEE-14BARA.txt           # Data and parameters for the IEEE 14-Bus system

├── IEEE-14BARA_SIMULATION.py # Simulation setup using Pandapower

└── README.md                 # Project documentation

## Prerequisites

- Python 3.8 or higher
- Pandapower library

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/emirhanersoy/CSA-ORPD.git
   cd CSA
2. Install dependencies:
   ```bash
   pip install pandapower

## Running the Project


1. Run the CSA algorithm:
   ```bash
   python CSA.py

2. (Optional)Simulate the IEEE 14-Bus System:
   ```bash
   python IEEE-14BARA_SIMULATION.py


# Features
Implementation of the Crystal Structure Algorithm for ORPD optimization
IEEE 14-Bus system simulation using Pandapower
Comprehensive constraint handling:

• Voltage levels
• Transformer tap settings
• Reactive power limits
