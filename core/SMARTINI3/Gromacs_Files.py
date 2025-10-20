#! /usr/bin/ipython3
import os
import numpy as np
import pandas as pd
import math
from   itertools import islice
import sys
import itertools 
from itertools import islice
import matplotlib.pyplot as plt



def MINIMIZATION(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    filename = f"Min_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.mdp"
    
    energy_groups = f"T-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                    f"H-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE}"
    
    energygrp_table = f"T-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                      f"T-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                      f"H-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                      f"H-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE}"
    
    with open(filename, "w") as file:
        file.write(
            "integrator               = steep\n"
            "emstep                   = 0.01\n"
            "emtol                    = 10\n"
            "nsteps                   = 1000000\n"
            "nstenergy                = 500\n"
            "nstlog                   = 500\n"
            "nstxout-compressed       = 1000\n"
            "; Single-range cutoff scheme\n"
            "cutoff-scheme            = group\n"
            "; PME electrostatics parameters\n"
            "coulombtype              = user        ; charges are set to zero in the topology file\n"
            f"rcoulomb                 = {rc}\n"
            "; Non-bonded parameters\n"
            "vdwtype                  = user\n"
            f"rlist                    = {rc}\n"
            f"rvdw                     = {rc}    ; LJ cut-off\n"
            f"energygrps               = {energy_groups}\n"
            f"energygrp-table          = {energygrp_table}\n"
        )




def NVT(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    filename = f"NVT_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.mdp"
    
    energy_groups = f"T-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                    f"H-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE}"
    
    energygrp_table = f"T-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                      f"T-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                      f"H-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE} " \
                      f"H-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE}"
    
    tc_grps = energy_groups
    
    with open(filename, "w") as file:
        file.write(
            "integrator               = sd         ; Langevin dynamics\n"
            "bd-fric                  = 0.5\n"
            "ld-seed                  = -1\n"
            "dt                       = 0.15\n"
            "nsteps                   = 20000       ; 15ns\n"
            "nstcomm                  = 1           ; Remove COM translational velocity every step\n"
            "; Output parameters\n"
            ";nstxout                  = 1\n"
            ";nstvout                  = 1\n"
            ";nstenergy                = 1\n"
            ";nstlog                   = 1\n"
            "; Single-range cutoff scheme\n"
            "cutoff-scheme            = group\n"
            "; PME electrostatics parameters\n"
            "coulombtype              = user        ; charges set to zero in topology\n"
            f"rcoulomb                 = {rc}\n"
            "couple-intramol           = yes\n"
            "; Non-bonded parameters\n"
            "vdwtype                  = user\n"
            f"rlist                     = {rc}\n"
            f"rvdw                      = {rc}    ; LJ cut-off\n"
            f"energygrps               = {energy_groups}\n"
            f"energygrp-table          = {energygrp_table}\n"
            "; Temperature coupling is on\n"
            "Tcoupl                   = no          ; handled by sd implicitly\n"
            f"tc_grps                  = {tc_grps}\n"
            "tau_t                     = 1.0    1.0\n"
            "ref_t                     = 315    315\n"
            "; Periodic boundary conditions\n"
            "pbc                       = xyz\n"
        )



    
    

