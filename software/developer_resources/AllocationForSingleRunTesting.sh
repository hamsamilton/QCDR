#!/bin/zsh

srun -N 1 -n 10 -A p31279 --mem=50G --partition=normal --time=08:00:00 --pty bash -l

