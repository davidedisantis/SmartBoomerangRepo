# Project Title

Dynamics and Control of Planetary Smart Boomerangs

## Description

The project focuses on the study of dynamical characteristics of boomerang and controls of Smart Boomerangs.
The work is divided into two folders: Standard Boomerangs and Smart Boomerangs.

The Standard Boomerangs folder contains three Excel sheets containing the aerodynamic coefficients of a NACA 0012 at 50.000 Re, one .m MATLAB file, and two Simulink models, one to run simulations without sensors and one to run simulations with an IMU sensor. 
The Boomerang_Sim file is the one you should open. There, dedicated sections allow you to model the boomerang and define its launch conditions.
Running the file allows you to choose the atmosphere where you want the boomerang to fly and which Simulink file to run.

The Smart Boomerangs folder contains three Excel sheets containing the aerodynamic coefficients of a NACA 0012 at 50.000 Re, one .m MATLAB file, and four Simulink models, one for each control method tested.
This codes are an expansion of the Standard Boomerangs one and shall be operated in a similar way. Running the SmartBoomerang_Sim.m code, lets you choose between the control methods available. In the Desired Trajectory section, a pre-determined path that you want the boomerang to follow shall be given as an input.
Dedicated sections within the code can be modified to tune the PID gains.

The default scripts refer to a 4-winged boomerangs, with launch characteristics that make it to circle back to the thrower in an undisturbed environment.
Each .m presents dedicated sections for plotting the results.
Some .mat files are available which contain different trajectories.
Dedicated sections allow to plot desired results.

## Getting Started

### Dependencies

The code was written in MATLAB/Simulink 2023a. Have an up-to-date version of Simulink to run the files.

### Executing program

Operate the .m files for running trajectory predictions and dynamical analyses.

## Help

For any issues, feel free to reach out to davide.disantis@gmail.com

## Author

Davide Di Santis
E-Mail: davide.disantis@gmail.com
LinkedIn: www.linkedin.com/in/davide-di-santis

## Acknowledgments
Thanks to Drs. Adrian Stoica and Marco B. Quadrelli for the supervision.
References for the Boomerang Dynamics are found in Flight Dynamics of the Boomerang, Part 1: Fundamental Analysis by Azuma et al.