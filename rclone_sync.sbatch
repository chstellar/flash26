#!/bin/bash
#
#SBATCH --job-name=rclone_sync
#SBATCH -p horence,normal
#SBATCH -c 4
#SBATCH --time=3:00:00

ml system
ml rclone
rclone copy results/ gdrive:Salzman\ Lab/Shared\ Materials/MetaSPLASH_all_results/