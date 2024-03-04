using Random
using Statistics
using Plots 
using CSV
using DataFrames

signal_acq = CSV.read("\\signal_acq.csv", DataFrame, header=false)
sig = CSV.read("\\sig.csv", DataFrame, header=false)


# Calcolo del Residual Sum of Squares (RSS)
rss = sum((Matrix(signal_acq) .- Matrix(sig)).^2);

# Calcolo del Total Sum of Squares (TSS)
tss = sum((Matrix(signal_acq) .- mean(Matrix(sig))).^2);

# Calcolo del coefficiente di determinazione R^2
r2 = (1 - rss / tss)*100
