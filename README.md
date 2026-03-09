# Theoretical Study of Stratified Photon Mapping with Kernel

## Overview

This repository contains my bachelor thesis titled **"Theoretical Study of Stratified Photon Mapping with Kernel"**, completed at Ludwig Maximilian University of Munich (LMU).

The work focuses on the theoretical analysis of stratified photon mapping within the framework of global illumination in physically based rendering. In contrast to many computer graphics projects that emphasize implementation, this thesis primarily investigates the mathematical and probabilistic structure underlying photon mapping algorithms.

## Abstract

In this work, the over- and under-estimation biases in the irradiance estimation of Photon Mapping, first reported in [Gar12] by Dr. García, are studied. In particular, these biases are analyzed for stratified Photon Mapping with kernels, adopting some of the stratification strategies used in [Har20].

For this purpose, the existence and behavior of these biases are examined for four kernels (Epanechnikov, Cone, Silverman, and Gaussian) used in [Gar14] and [Gan18]. It is mathematically and empirically demonstrated that the over- and under-estimation biases in stratified Photon Mapping cannot in general be eliminated using kernels; however, a reduction of these biases is possible when kernels are employed.

Additionally, the improvement in the standard deviation for the different kernels is empirically studied under various stratification strategies. A stratification of the scene is simulated by stratifying the light source, aiming to generate a secondary stratification over the scene. In particular, both the position and the emission direction of the photons are stratified.

## Topics

- Computer Graphics  
- Global Illumination  
- Photon Mapping  
- Monte Carlo Rendering  
- Light Transport Simulation  
- Kernel Density Estimation  
  

## Acknowledgments

I would like to thank **Dr. Rubén Jesús García-Hernández** for his supervision and for providing the reference code implementation used in the experimental part of this work.