#!/usr/bin/env bash

# ---------------------------
# Script for post-processing GROMACS simulations
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

# ---------------------------
# File name templates
# ---------------------------
npt_trr="npt_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.trr"
npt_tpr="npt_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.tpr"
index_file="index-${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.ndx"
conf_file="conf_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}_.gro"

# ---------------------------
# Extract snapshots and fix PBC
# ---------------------------
echo "0" | gmx trjconv -s "$npt_tpr" -f "$npt_trr" -o "$conf_file" -sep -pbc whole

# ---------------------------
# Box dimensions for Area Compressibility
# ---------------------------
for dim in LX LY LZ; do
    idx=$([[ $dim == "LX" ]] && echo 12 || [[ $dim == "LY" ]] && echo 13 || echo 14)
    gmx energy -f "${npt_trr%.trr}.edr" -o "${dim}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" -b 0 -e 252000 <<< "$idx"
done

# ---------------------------
# Density profile
# ---------------------------
gmx density -f "$npt_trr" -n "$index_file" -s "$npt_tpr" \
    -o "density_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" \
    -ng 3 -sl 1000 <<< "0 2 3"

# Normalize density
tail -n +27 "density_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" > "New_density_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg"
awk '{print $1/8.69467, "\t", $4/774.043}' "New_density_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" > "Normalised_density_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg"

# ---------------------------
# Radial Distribution Function (RDF)
# ---------------------------
rdf_types=( "T_T" "H_H" "T_H" )
refs=(2 3 2)
sels=(2 3 3)

for i in ${!rdf_types[@]}; do
    gmx rdf -f "$npt_trr" -n "$index_file" -s "$npt_tpr" \
        -o "g_${rdf_types[i]}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" \
        -seltype mol_com -bin 0.01 -ref ${refs[i]} -sel ${sels[i]}
    
    # Normalize RDF
    tail -n +27 "g_${rdf_types[i]}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" > "New_g_${rdf_types[i]}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg"
    awk '{print $1/4.070, "\t", $2/2.7}' "New_g_${rdf_types[i]}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" > "Normalised_g_${rdf_types[i]}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg"
done

# ---------------------------
# Line tension calculation (pressure components)
# ---------------------------
p_components=(21 25 29)
dims=(Pxx Pyy Pzz)
for i in ${!dims[@]}; do
    gmx energy -f "nvt_LT_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.edr" \
        -o "${dims[i]}_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.xvg" \
        -b 0 -e 252000 <<< "${p_components[i]}"
done

# ---------------------------
# Compute metrics using Python scripts
# ---------------------------
results_file="RESULTS_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}.txt"
> "$results_file"

for script in AC DP LT OPT_DP OPT_RDF BendingM; do
    output=$(python "./${script}.py" $epsilon $rmin $r1 $rc $bFC $aFC $Friction $NFE)
    echo "$output" >> "$results_file"
done

# ---------------------------
# Clean up
# ---------------------------
mv "conf_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}_999.gro" \
   "CONF_${epsilon}_${rmin}_${r1}_${rc}_${bFC}_${aFC}_${Friction}_${NFE}_999.gro"
rm conf*.gro
rm *.cpt *.pdb New_g_T_T_*.xvg Normalised_density_*.xvg New_density_*.xvg Normalised_g_T_T_*.xvg NPT_PP_*.xvg
rm \#*

# ---------------------------
# Time calculation
# ---------------------------
time_file="TIME.Threetxt${NFE}"
> "$time_file"
elapsed=$(($SECONDS - $start_time))
echo "$(format_time $elapsed)" >> "$time_file"

# ---------------------------
# File management: move results
# ---------------------------
Main_Dir="/home/uni08/soleimani/RUNS/Thesis/Membrane_Project/PAPER/9/"
mkdir -p "$Main_Dir/FIGURES" "$Main_Dir/EDR" "$Main_Dir/TIME"

cp *.png "$Main_Dir/FIGURES"
cp *.edr "$Main_Dir/EDR"
cp TIME.Onetxt${NFE} TIME.Twotxt${NFE} TIME.Threetxt${NFE} "$Main_Dir/TIME"
