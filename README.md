# Conceptual Aircraft Design Optimization Solver (CADOS)
Course project for ENGR493 Intro to Aircraft Design/AIAA Carrier-based Strike Fighter Competition Project

Welcome to the project!
The general structure is as follows:

## .m Files
1. testerexpanded.m ~ This file is the main looping function
2. aa.m ~ This file contains all aerodynamic relationships
3. atmosphere.m ~ This file contains useful atmospheric parameters
4. const.m ~ This file should be used for constants, i.e. g=9.81
5. performance.m ~ This file contains all performance relationships
6. subsystems.m ~ This file contains all subsystem calculations
7. weight.m ~ This file contains aircraft total & fuel weight estimations

## .mat/.csv Files
1. /data/* ~ These files contain airfoil data, formatted as xf-airfoiltype-il-Re#.csv
2. a.mat ~ This file contains data for supersonic drag calculations
3. airfoils.mat ~ These files contain the airfoil data from /data/ but in .mat form
4. Atmosphere.csv ~ This file contains relevant atmospheric data
5. output.mat ~ This is the output file used for the original report

## Others
1. .git ~ This folder is related to github stuff, don't worry about it
2. .gitignore ~ same as above