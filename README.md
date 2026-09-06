# brownian-particle-tracking-analysis
Pipeline for 2D Brownian particle tracking and trajectory analysis.

This repository contains the MATLAB code I worked on during my research internship at IMFT (Institut de Mécanique des Fluides de Toulouse) / CNRS. 

The goal of this project was to build a full pipeline to process raw, high-frequency experimental data, filter out the noise, and extract statistical metrics from particle trajectories.

# Repository Structure

Here is a quick explanation of how the code is organized:

* **/pipeline**: Contains the main scripts to run the tracking process from start to finish (including `submission_fulltracking.m`). This is where the raw data ingestion and initial filtering happen.
* **/analysis**: My core statistical scripts. This includes `analyzeTrajectory.m` and `vitesse.m`, which I wrote to compute Mean Squared Displacement (MSD), velocity distributions, and generate PDFs.
* **/DC-MSS-JAQAMAN-VEGA-UT-SOUTHWESTERN**: An algorithm for motion regime detection originally developed by UT Southwestern, which I integrated and adapted for my specific datasets.
* **/external**: Other third-party dependencies and helper functions used across the project.
