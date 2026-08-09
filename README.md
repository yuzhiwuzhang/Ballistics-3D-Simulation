# Ballistics-3D-Simulation
3D artillery ballistic simulation with real-time trajectory visualization.  Includes aerodynamic drag (Mach-dependent), wind field, Coriolis effect, Magnus effect,  and standard atmosphere model. Developed in Python + C for embedded deployment.  TRL-3 proof-of-concept for future fire control research.  中文简述：含气动阻力、风场、科里奥利力、马格努斯效应的3D弹道仿真系统。
# Ballistics-3D-Simulation

3D artillery ballistic simulation with real-time trajectory visualization.

## Features
- 3D trajectory with Mach-dependent drag coefficient
- International Standard Atmosphere (ISA) model
- Wind field + Coriolis effect + Magnus effect
- Comparative analysis (with/without lateral forces)
- Python visualization + C version for embedded deployment
- TRL-3 concept validation

## Project Structure
- `3D弹道模型.py` — Main simulation + visualization
- `ballistics.c` — Embedded C version for real-time deployment

## Results
- 3D trajectory plot with rotatable view
- Top-down view showing lateral deviation
- Side view for range and altitude
- Time-history of lateral offset

## Physics Model
- Mach-dependent drag coefficient (0–3 Mach)
- Altitude-dependent density/temperature/speed of sound
- Gravity variation with altitude
- Coriolis force
- Magnus effect
- 3D wind field

## Future Work
- RK4 integrator
- 6DOF rigid body model
- Real-time sensor data interface

## Author
Independent researcher, high school senior (Class of 2027)
