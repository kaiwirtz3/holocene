# master bash script for the joint analysis of archaeological C14 and paleo-climate proxy data
# kai wirtz Dec 2023
# edit scdir='out/'; and 'addpath('~/tools/m_map');' in load_pars.m

#  reads C14_europe_*  saves eurogrid_*.bin
Rscript grid_growth.r
# on HPC: sbatch --array=1-64 slurm.sh

Rscript collect.r # loads bin/eurogrid_1-64, writes bin/eurodat_all
Rscript cluster.r # in parallel mode slumclust.sh
#  forms regions using an extended kmeans algorithm with optimized total number of clusters (for each RegionPatch)
#  saves clusti.mat

# region kriging on a grid  - based on cluster points
matlab -nodesktop -r "try; make_grid; catch; end; quit" > LogFile2 2> ErrorLogFile2

Rscript spd_growth.r  # reads clusti, writes AllPop for time slices
#or sbatch --array=0-101 slurmspd.sh (for 6 SPD settings: 6*17=102, start with 0!)

# pooled method
Rscript spd_pooled.r #reads C14_europe0 (or C14_EA, C14_NIreland), writes AllPop_all

# collects RGR from time slices
plot_ContRGR %reads AllPop_i writes AllPop_tag_all avg_rgr_ (avg_rgr_all all major methods)  RGR_MethComp.png

#plot_varmap_slice %reads AllPop_i

# -------------------------------------
# process climate proxy data
proxy_dtw  # applies DTW writes dtw_proxydata
## pca_proxydata #
#? prepare_climdata # reads events/proxydescription_488, dtw_proxydata saves 'proxydat_*'

collect_ts # integrate and smooth most time-series

glmloop #  runs GLM for boom/bust probability with 1-4 input variables

# prepare time series
check_stat # reads avg_rgr_all archoccdens Us16_comp bog_std dtwpca/dtwpca_proxydata_; writes  target_ts_
           # creates pca-based index 'climate stability', smoothes RGR with movweighavg(...,202,50);
# combine climate and RGR
#?? show_stat # reads target_ts_67 ; writes show_stat_.png glm_coeff_.tex

plot_stat # reads target_ts_13 ; writes show_stat_.png glm_coeff_.tex
