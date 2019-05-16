#!/bin/bash -l

#SBATCH --job-name=splashback_%j
#SBATCH --time=03:00:00
##SBATCH --ntasks=14
#SBATCH --exclusive
#SBATCH --partition=broadwl
#SBATCH --account=pi-chihway
#SBATCH --output=log/splashback-%A_%a.out
#SBATCH --error=log/splashback-%A_%a.err
#SBATCH --array=0-99%100

data_dir=/project/kicp/chihway/brutus/splashback/

source /project2/chihway/setup/setup_midway2.sh

jkid=${SLURM_ARRAY_TASK_ID}

python measure_Sigmag_sdss.py 0.1 0.2 4 -100 -19.43 20 100 ${data_dir}data/sdss/sdss_redmapper_jk_allz.fits ${data_dir}data/sdss/sdss_redmapper_randoms_jk_allz.fits ${data_dir}data/sdss/sdss_gal_morph_jk.fits ${data_dir}data/sdss/sdss_gal_ran_small_jk.fits $jkid 9930.94 results/l20_100_sdss/Sigmag_${jkid}.npz 25


