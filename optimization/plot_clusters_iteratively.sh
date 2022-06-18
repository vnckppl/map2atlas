#!/bin/bash
# * Usage
usage () {
    cat<<EOF
    Mandatory options:
    -f path to map2atlas.py script
    -o output folder where data will be written to
    -m map (nifti file) of area you want to grab clusters from your atlas for
    -a altas (nifti file) with integers representing areas of the same cluster
    -l threshold
    -h threshold
    -p indicate if you to run both threshold1 and threshold2

    Optional options:
    -s take a screenshot of the selected ROIs
EOF
}

# * Input
while getopts "f:o:m:a:l:h:p:s" OPTION
do
    case $OPTION in
        f)
            map2atlas_function="${OPTARG}"
            ;;
        o)
            odir="${OPTARG}"
            ;;
        m)
            map="${OPTARG}"
            ;;
        a)
            atlas="${OPTARG}"
            ;;
        l)
            mthrl="${OPTARG}"
            ;;
        h)
            mthrh="${OPTARG}"
            ;;
        p)
            parts="${OPTARG}"
            ;;
        s)
            screenshot=1
            ;;
        ?)
        usage
        exit
        ;;
    esac
done


# * Test input
if [ -z "${1}" ]; then
    usage
fi
if [ -z ${screenshot} ]; then
    screenshot=0
fi

# * Function for calculating proportion of overlap
# Proportion overlap between a) selected voxels inside the selection mask and
# total voxels inside the mask; and b) selected voxels outside the selection
# mask and total voxels outside the mask.
in_out_calc() {

    # * Environment
    base="${1}"
    tdir=${base}/tmp
    mkdir -p "${tdir}"

    # ** Input files
    # Neurosynth mask (bin)
    mask=${base}/map_t.nii.gz
    # Full Shen atlas (all regions)
    atlas=${base}/atlas_r.nii.gz
    # Selected clusters (all regions)
    clusters=${base}/ROI_selection.nii.gz

    # * Binarize input files
    atlas_bin=${tdir}/atlas_r_bin.nii.gz
    fslmaths "${atlas}" -bin "${atlas_bin}"

    cl_bin=${tdir}/ROI_selection_bin.nii.gz
    fslmaths "${clusters}" -bin "${cl_bin}"

    # * Create inverse of the neurosynth mask
    mask_inv=${tdir}/map_t_inv.nii.gz
    fslmaths "${mask}" -binv "${mask_inv}"

    # * Create maps for areas of the volume and grab the voxels
    # 1) Voxels of all AVAILABLE atlas parcels INSIDE the neurosynth mask
    a_in_n=${tdir}/a_in_n.nii.gz
    fslmaths "${atlas_bin}" -mas "${mask}" "${a_in_n}"
    a_in_n_vox=$(fslstats "${a_in_n}" -V | awk '{ print $1 }')

    # 2) Voxels of all AVAILABLE atlas parcels OUTSIDE the neurosynth mask
    a_out_n=${tdir}/a_out_n.nii.gz
    fslmaths "${atlas_bin}" -mas "${mask_inv}" "${a_out_n}"
    a_out_n_vox=$(fslstats "${a_out_n}" -V | awk '{ print $1 }')

    # 3) Voxels of all SELECTED atlas parcels INSIDE the neurosynth mask
    s_in_n=${tdir}/s_in_n.nii.gz
    fslmaths "${cl_bin}" -mas "${mask}" "${s_in_n}"
    s_in_n_vox=$(fslstats "${s_in_n}" -V | awk '{ print $1 }')

    # 4) Voxels of all SELECTED atlas parcels OUTSIDE the neurosynth mask
    s_out_n=${tdir}/s_out_n.nii.gz
    fslmaths "${cl_bin}" -mas "${mask_inv}" "${s_out_n}"
    s_out_n_vox=$(fslstats "${s_out_n}" -V | awk '{ print $1 }')

    # * Calculate the ratios
    perc_inside=$(echo "scale=4; ${s_in_n_vox} / ${a_in_n_vox}" | bc -l)
    perc_outside=$(echo "scale=4; ${s_out_n_vox} / ${a_out_n_vox}" | bc -l)
    diff=$(echo "scale=4; ${perc_inside} - ${perc_outside}" | bc -l)

    # * Export
    ofile=${base}/optimization.csv
    echo "p_in,p_out,diff" > "${ofile}"
    echo "${perc_inside},${perc_outside},${diff}" >> "${ofile}"

    # * Echo results
    echo "${perc_inside},${perc_outside},${diff}"
}


# * Environment
mkdir -p "${odir}"

## ** Output file for counting number of selected clusters per iteration
ofile1=${odir}/counts.csv
echo "othr1,othr2,cluster_count" > "${ofile1}"

## ** Output file for optimization values
ofile2=${odir}/optimization.csv
echo "othr1,othr2,prop_inside,prop_outside,difference" > "${ofile2}"


## * Iterate over values of othr1
time (
    for othr1 in $(seq -w 5 5 100); do

          ## ** Iterate over values of othr2
          for othr2 in $(seq -w 5 5 100); do

              ## *** Announce
              echo "Threshold 1 iteration: ${othr1}/100 -- Threshold 2 iteration: ${othr2}/100"
              
              ## *** Output folder
              odir2=${odir}/${othr2}/${othr1}
              mkdir -p "${odir2}"

              ## *** Run map2atlas    
              python3 "${map2atlas_function}" \
                      "${map}" \
                      "${atlas}" \
                      "${odir2}" \
                      --othr1 "${othr1}" \
                      --othr2 "${othr2}" \
                      --mthr "${mthrl}" "${mthrh}" \
                      --parts "${parts}" \
                      --com \
                  &> "${odir2}/log.txt"

              ## *** Count clusters
              n=$(sed 's/,/ /g' "${odir2}/ROI_selection_list.csv" | wc -w)
              echo "${othr1},${othr2},${n}" >> "${ofile1}"

              ## *** Run calculation of voxels inside and outside the mask
              ioc=$(in_out_calc "${odir2}")
              echo "${othr1},${othr2},${ioc}" >> "${ofile2}"

              ## *** Create screenshot
              if [ "${screenshot}" -eq 1 ]; then
                  fsleyes \
                      render \
                      -of "${odir2}/ss.png" \
                      --size 1280 800 \
                      -ds world \
                      "${FSLDIR}"/data/standard/MNI152_T1_0.5mm.nii.gz \
                      "${odir2}"/atlas_r.nii.gz \
                          -ot label --lut random_big -w 0 -n Shen256 -a 25 \
                      "${odir2}"/ROI_selection.nii.gz \
                          -ot label --lut random_big -w 0 -n Selected_ROIs \
                      "${odir2}"/map_t \
                          -ot label -o -w 3 --lut random \
                      &>> "${odir2}/log.txt"
              fi
          done
    done
)

# * Run R script to plot output and find optimal combination of othr1 and othr2
if [ "$(command -v Rscript)" ]; then
    sdir=$(dirname "${0}")
    Rscript "${sdir}/plot_clusters_iteratively.R" "${odir}"
else
    echo "Rscript not found. Make sure you have R with the ggplot2 package installed."
    echo "If R is installed, make sure that the executive is added to your PATH variable."
    exit 1
fi
