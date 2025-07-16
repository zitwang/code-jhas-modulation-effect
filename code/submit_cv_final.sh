#!/usr/bin/env bash
#
# Usage: ./submit_cv_final.sh <n_cv_folds> [mri_role] [MRI_effect]
#    e.g. ./submit_cv_final.sh 5 both binary

n_cv_folds=${1:?Please supply number of folds}
mri_role=${2:-both}
MRI_effect=${3:-binary}

# adjust these to taste
workdir="/users/zwang3/PAS_final"
partition="shared"
mem="5G"
cpus=1
email="zwang238@jhmi.edu"

for ((to_mask=1; to_mask<=n_cv_folds; to_mask++)); do
  job_name="cv${to_mask}_${n_cv_folds}_${mri_role}_${MRI_effect}"

  sbatch \
    --partition="$partition" \
    --mem="$mem" \
    --cpus-per-task="$cpus" \
    --job-name="$job_name" \
    --output="${workdir}/logs/${job_name}_%j.out" \
    --error="${workdir}/logs/${job_name}_%j.err" \
    --mail-user="$email" \
    --mail-type=END,FAIL \
    --export=K="$n_cv_folds",to_mask="$to_mask",mri_role="$mri_role",workdir="$workdir",MRI_effect="$MRI_effect",job_name="$job_name" \
    submit_single_cv_final.sh
done
