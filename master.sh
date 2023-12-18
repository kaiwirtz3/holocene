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
##make_grid.m

Rscript spd_growth.r  # reads clusti, writes AllPop for time slices
#or sbatch --array=0-101 slurmspd.sh (for 6 SPD settings: 6*17=102, start with 0!)

# pooled method
Rscript spd_all.r #reads C14_europe_neo, writes AllPop_all

# collects RGR from time slices
plot_ContRGR %reads AllPop_i writes AllPop_tag_all avg_rgr_ (avg_rgr_all all major methods)  RGR_MethComp.png

# plot map with proxy locations and occupation boxes
plot_sites

plot_varmap_slice %reads AllPop_i

# plots and saves regional RGR ; e.g., writes Ire_rgr
extract_rgr

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

#cmp_plot_Growth # loads avg_rgr_all or archoccdens (from occdens)
                # with pvi=[10 11] compares to regional boom/busts

# compares with SPD RGR  South America from Riris et al 2019
#plot_spd_SA
#write_ts

matlab -nodesktop -r "try; master2; catch; end; quit" > LogFile2 2> ErrorLogFile2
