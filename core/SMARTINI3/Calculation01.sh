#!/usr/bin/env bash

# ---------------------------
# Script to run GROMACS simulations
# ---------------------------

# Start timer
start_time=$SECONDS

# Read input arguments
epsilon=$1
rmin=$2
r1=$3
rc=$4
bFC=$5
aFC=$6
Friction=$7
NFE=$8

# Number of cores
Core_NUMBER=4

# Directory output file
dir_file="Directory_${Core_NUMBER}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.txt"
> "$dir_file"
echo "$(pwd)" >> "$dir_file"

# Function to format time as HH:MM:SS
format_time() {
  local total_seconds=$1
  ((h=total_seconds/3600))
  ((m=(total_seconds%3600)/60))
  ((s=total_seconds%60))
  printf "%02d:%02d:%02d\n" $h $m $s
}

# Index file name
index_file="index-${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.ndx"

# ---------------------------
# Create index groups
# ---------------------------
gmx make_ndx -f Bilayer.gro -o "$index_file" << EOF
name 2 T-${epsilon}-${rmin}-${r1}-${rc}-${bFC}-${aFC}-${Friction}-${NFE}
name 3 H-${epsilon}-${rmin}-${r1}-${rc}-${bFC}-${aFC}-${Friction}-${NFE}
q
EOF

# ---------------------------
# Energy Minimization
# ---------------------------
em_mdp="Min_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.mdp"
em_tpr="em_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.tpr"
em_gro="em_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.gro"

gmx grompp -f "$em_mdp" -c Bilayer.gro -n "$index_file" \
    -p "topol_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.top" \
    -o "$em_tpr" -maxwarn -1

gmx mdrun -s "$em_tpr" -table table.xvg -o "$em_gro" -e "$em_gro" -x "$em_gro" \
    -c "$em_gro" -cpo "$em_gro" -nt $Core_NUMBER

# ---------------------------
# NVT Equilibration
# ---------------------------
nvt_mdp="NVT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.mdp"
nvt_tpr="nvt_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.tpr"
nvt_gro="nvt_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.gro"

gmx grompp -f "$nvt_mdp" -c "$em_gro" -n "$index_file" \
    -p "topol_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.top" \
    -o "$nvt_tpr" -maxwarn -1

gmx mdrun -s "$nvt_tpr" -table table.xvg -o "$nvt_gro" -e "$nvt_gro" -x "$nvt_gro" \
    -c "$nvt_gro" -cpo "$nvt_gro" -nt $Core_NUMBER

# ---------------------------
# NPT Equilibration
# ---------------------------
npt_mdp="NPT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.mdp"
npt_tpr="npt_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.tpr"
npt_gro="npt_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.gro"

gmx grompp -f "$npt_mdp" -c "$nvt_gro" -n "$index_file" \
    -p "topol_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.top" \
    -o "$npt_tpr" -maxwarn -1

gmx mdrun -s "$npt_tpr" -table table.xvg -o "$npt_gro" -e "$npt_gro" -x "$npt_gro" \
    -c "$npt_gro" -cpo "$npt_gro" -nt $Core_NUMBER

# ---------------------------
# Time calculation
# ---------------------------
time_file="TIME.Onetxt${NFE}"
> "$time_file"
elapsed=$(($SECONDS - $start_time))
echo "$(format_time $elapsed)" >> "$time_file"
