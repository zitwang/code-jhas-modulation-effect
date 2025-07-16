#!/usr/bin/env bash
#
# ── SLURM directives ───────────────────────────────────────────────────────────
#SBATCH --time=48:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zwang238@jhmi.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL            # imports K,to_mask,mri_role,workdir,job_name

# ── payload ────────────────────────────────────────────────────────────────────

# sanity-check that the five variables came through
: "${workdir:?need workdir}"
: "${job_name:?need job_name}"
: "${K:?need K}"
: "${to_mask:?need to_mask}"
: "${mri_role:?need mri_role}"
: "${MRI_effect:?need MRI_effect}"

# ensure logs directory exists
mkdir -p "${workdir}/logs"

echo "[$(date)] Starting ${job_name} on $(hostname)"

module load conda_R

# run your R script; capture everything into its own log
Rscript "${workdir}/code/cv-to-run-single.R" \
  "${K}" "${to_mask}" "${mri_role}" "${workdir}" "${MRI_effect}"\
  &> "${workdir}/logs/${job_name}_rscript_log.txt"
