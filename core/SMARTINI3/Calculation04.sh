#!/usr/bin/bash

# Track elapsed time
secs=$SECONDS

# Function to format seconds into HH:MM:SS
format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
}

# ------------------------------
# Script arguments
# $1: epsilon
# $2: rmin
# $3: r1
# $4: rc
# $5: bFC
# $6: aFC
# $7: FRIC
# $8: NFE
# $9: Temp
# ------------------------------

CORE_NUMBER=4

# ------------------------------
# Directory log
# ------------------------------
DIR_LOG="Directory_${CORE_NUMBER}_$1_$2_$3_$4_$5_$6_$7_$8_$9.txt"
> "$DIR_LOG"
echo "$(pwd)" >> "$DIR_LOG"

# ------------------------------
# Index file
# ------------------------------
INDEX_FILE="index-$1_$2_$3_$4_$5_$6_$7_$8.ndx"
gmx make_ndx -f system_mix.gro -o "$INDEX_FILE" << EOF
name 2  T-$1-$2-$3-$4-$5-$6-$7-$8
name 3  H-$1-$2-$3-$4-$5-$6-$7-$8
q
EOF

# ------------------------------
# Filenames (for cleaner code)
# ------------------------------
MIN_MDP="Min_$1_$2_$3_$4_$5_$6_$7_$8_$9.mdp"
NVT_MDP="NVT_$1_$2_$3_$4_$5_$6_$7_$8_$9.mdp"
NPT_MDP="NPT_$1_$2_$3_$4_$5_$6_$7_$8_$9.mdp"

TOPOL="topol_$1_$2_$3_$4_$5_$6_$7_$8_$9.top"

EM_TPR="em_$1_$2_$3_$4_$5_$6_$7_$8.tpr"
NVT_TPR="nvt_$1_$2_$3_$4_$5_$6_$7_$8.tpr"
NPT_TPR="npt_$1_$2_$3_$4_$5_$6_$7_$8.tpr"

# ------------------------------
# Energy minimization
# ------------------------------
gmx grompp -f "$MIN_MDP" -c system_mix.gro -n "$INDEX_FILE" -p "$TOPOL" -o "$EM_TPR" -maxwarn -1
gmx mdrun -s "$EM_TPR" -table table.xvg -o "${EM_TPR%.tpr}" -e "${EM_TPR%.tpr}" -x "${EM_TPR%.tpr}" \
          -c "${EM_TPR%.tpr}" -cpo "${EM_TPR%.tpr}" -nt "$CORE_NUMBER"

# ------------------------------
# NVT equilibration
# ------------------------------
gmx grompp -f "$NVT_MDP" -c "${EM_TPR%.tpr}.gro" -n "$INDEX_FILE" -p "$TOPOL" -o "$NVT_TPR" -maxwarn -1
gmx mdrun -s "$NVT_TPR" -table table.xvg -o "${NVT_TPR%.tpr}" -e "${NVT_TPR%.tpr}" -x "${NVT_TPR%.tpr}" \
          -c "${NVT_TPR%.tpr}" -cpo "${NVT_TPR%.tpr}" -nt "$CORE_NUMBER"

# ------------------------------
# NPT equilibration
# ------------------------------
gmx grompp -f "$NPT_MDP" -c "${NVT_TPR%.tpr}.gro" -n "$INDEX_FILE" -p "$TOPOL" -o "$NPT_TPR" -maxwarn -1
gmx mdrun -s "$NPT_TPR" -table table.xvg -o "${NPT_TPR%.tpr}" -e "${NPT_TPR%.tpr}" -x "${NPT_TPR%.tpr}" \
          -c "${NPT_TPR%.tpr}" -cpo "${NPT_TPR%.tpr}" -nt "$CORE_NUMBER"

# ------------------------------
# Time tracking
# ------------------------------
TIME_LOG="TIME.FIVEtxt$9"
> "$TIME_LOG"
echo "$(format_time $SECONDS)" >> "$TIME_LOG"
