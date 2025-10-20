#!/usr/bin/env bash

# ---------------------------
# Script to run GROMACS Line Tension simulations
# ---------------------------

# Start timer
start_time=$SECONDS

# ---------------------------
# Function to format time as HH:MM:SS
# ---------------------------
format_time() {
  local total_seconds=$1
  ((h=total_seconds/3600))
  ((m=(total_seconds%3600)/60))
  ((s=total_seconds%60))
  printf "%02d:%02d:%02d\n" $h $m $s
}

# ---------------------------
# Input arguments
# ---------------------------
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

# ---------------------------
# Index file
# ---------------------------
index_file="index-LT-${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.ndx"

# Create index groups
gmx make_ndx -f Bilayer_LT_new.gro -o "$index_file" << EOF
name 2 T-${epsilon}-${rmin}-${r1}-${rc}-${bFC}-${aFC}-${Friction}-${NFE}
name 3 H-${epsilon}-${rmin}-${r1}-${rc}-${bFC}-${aFC}-${Friction}-${NFE}
q
EOF

# ---------------------------
# Energy Minimization
# ---------------------------
em_mdp="Min_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.mdp"
em_tpr="em_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.tpr"
em_gro="em_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.gro"

gmx grompp -f "$em_mdp" -c Bilayer_LT_new.gro -n "$index_file" \
    -p "topol_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.top" \
    -o "$em_tpr" -maxwarn -1

gmx mdrun -s "$em_tpr" -table table.xvg -o "$em_gro" -e "$em_gro" -x "$em_gro" \
    -c "$em_gro" -cpo "$em_gro" -nt $Core_NUMBER

# ---------------------------
# NVT Equilibration
# ---------------------------
nvt_mdp="NVT_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.mdp"
nvt_tpr="nvt_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.tpr"
nvt_gro="nvt_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.gro"

gmx grompp -f "$nvt_mdp" -c "$em_gro" -n "$index_file" \
    -p "topol_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.top" \
    -o "$nvt_tpr" -maxwarn -1

gmx mdrun -s "$nvt_tpr" -table table.xvg -o "$nvt_gro" -e "$nvt_gro" -x "$nvt_gro" \
    -c "$nvt_gro" -cpo "$nvt_gro" -nt $Core_NUMBER

# ---------------------------
# Time calculation
# ---------------------------
time_file="TIME.Twotxt${NFE}"
> "$time_file"
elapsed=$(($SECONDS - $start_time))
echo "$(format_time $elapsed)" >> "$time_file"
