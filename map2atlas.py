#! /usr/bin/env python3

# * map2atlas
# Vincent Koppelmans
# 2021-04-02

# * Libraries
import argparse
import os
import nibabel as nb
import numpy as np
from scipy import ndimage
from nipype.interfaces import freesurfer, fsl

# * Import arguments
# Gather arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='map2atlas')

    parser.add_argument('mapf',
                        help='The input image for which you want to know which '
                        'atlas ROIs it overlaps with.')

    parser.add_argument('atlas',
                        help='Atlas from which you want to know which ROIs '
                        'overlap with your map file.')

    parser.add_argument('outdir',
                        help='Output folder')

    parser.add_argument('--othr1',
                        help='Minimal percent overlap of the ROI with the map '
                        'to be selected. E.g., an ROI needs to be filled for '
                        'at least x percent by the map before it will be '
                        'selected. Default is 50 percent.',
                        default=50,
                        choices=range(0, 101),
                        metavar='[0-100]',
                        type=float)

    parser.add_argument('--othr2',
                        help='Minimal percent overlap of the map clusters and '
                        'the ROI to be selected. Default is 50 percent. E.g., '
                        'After the map is broken up into clusters, each '
                        'needs to fit for at least x percent inside an ROI '
                        'before that ROI is selected. Default = 50 percent',
                        default=50,
                        choices=range(0, 101),
                        metavar='[0-100]',
                        type=float)

    parser.add_argument('--mthr',
                        help='Map threshold (lower, upper) level to create a '
                        'binary image from your map file. Default = 0,1.',
                        nargs=2,
                        default=[0, 1],
                        metavar=('low', 'high'))

    parser.add_argument('--com',
                        help='Only select ROIs if their center-of-mass is '
                        'inside the map.',
                        default=False,
                        action='store_true')

    parser.add_argument('--parts',
                        help='Run either part1, or part 1 and part2 of this '
                        'script (see under description). Default = both',
                        default="both",
                        choices=["p1", "both"])

    args = parser.parse_args()


# * Empty Object Creator
class Scratch(object):
    pass


# * For testing
testing = False
if testing:
    args = Scratch()

    # Default mode network Zmap
    args.mapf = ('/Users/vincent/Data/Documents/Utah/Kladblok'
                 '/20170616_Weisenbach/20170822_Suicide'
                 '/20170823_Suicide_Proposal_Prelim/03_WFND/20_components'
                 '/Images_for_Proposal/2D/DMN.nii.gz')

    # Shen Atlas in MNI space
    args.atlas = ('/Users/vincent/Data/Documents/Utah/Kladblok'
                  '/20171002_Neuroimaging/20170629_Atlases'
                  '/20200107_Shen_in_MNI152/output/shenToMNI152.nii.gz')

    args.outdir = ('/Users/vincent/Data/tmp/20210402_map2atlas_test')

    args.othr1 = 20
    args.othr2 = 10
    args.mthr = [1, 100]
    args.com = True
    args.parts = "both"


# * Build map2atlas Object
class map2atlas:
    def __init__(self):

        # ** Environment
        class arguments:
            def __init__(self):
                self.mapf = args.mapf
                self.atlas = args.atlas
                self.outdir = args.outdir
                self.othr1 = float(args.othr1)
                self.othr2 = float(args.othr2)
                self.mthr = args.mthr
                self.com = args.com
                self.parts = args.parts
        # Call: Prepare map data
        self.args = arguments()

        # ** Create output folder
        os.makedirs(self.args.outdir, exist_ok=True)

        # ** Prepare map data
        class map_data:
            def __init__(self, mapf, mthr, outdir):

                # *** Load map data
                self.img = nb.load(mapf)

                # *** Threshold and binarize map data
                self.binimg = self.img.get_fdata()
                self.binimg = np.where(
                    self.binimg < float(mthr[1]), self.binimg, 0
                )
                self.binimg = np.where(self.binimg > float(mthr[0]), 1, 0)

                # *** Save thresholded image
                self.outimg = nb.Nifti1Image(
                    self.binimg, self.img.affine, self.img.header
                )
                self.outfile = outdir + '/map_t.nii.gz'
                nb.save(self.outimg, self.outfile)
        # Call: Prepare map data
        self.map_data = map_data(
            self.args.mapf, self.args.mthr, self.args.outdir
        )

        # ** Preapare atlas data
        class atlas_data:
            def __init__(self, odir, mapf, atlas):

                # *** Load atlas data
                self.img = nb.load(atlas)

                # *** Reslice atlas image to map space
                # nipype object for reslicing
                self.outfile = odir + '/atlas_r.nii.gz'
                self.reslice = freesurfer.MRIConvert()
                self.reslice.inputs.in_file = atlas
                self.reslice.inputs.out_file = self.outfile
                self.reslice.inputs.reslice_like = mapf
                self.reslice.inputs.resample_type = 'nearest'
                self.reslice.terminal_output = 'none'
                self.results = self.reslice.run()

                # *** List unique ROI values in the atlas file
                self.imgrdata = nb.load(self.outfile)
                self.roinums = np.unique(self.imgrdata.get_fdata())
                # Remove zero, because that is not a label
                self.roinums = self.roinums[self.roinums != 0]

        # Call: Prepare atlas data
        self.atlas_data = atlas_data(
            self.args.outdir, self.args.mapf, self.args.atlas
        )

        # ** -------------------------------------------------------------------
        # ** PART 1: Select ROIs that sufficiently fit inside the map image
        # ** -------------------------------------------------------------------
        if args.parts == "p1" or args.parts == "both":

            # ** Announce
            print("PART 1: Selecting ROIs that sufficiently fit inside the map"
                  "image...")

            # ** Function to test if atlas ROI overlaps with the map
            def test_roi(i, atld, mapd, thr):

                # *** Extract ROI (binary)
                roibin = np.where(atld == i, 1, 0)

                # *** Size of the ROI in voxels
                roisize = np.sum(roibin)

                # *** If the center of mass flag was enabled
                if self.args.com:

                    # Calculate center of mass coordinates
                    roicom = ndimage.center_of_mass(roibin)

                    # Get the value of the COM coordinate from the map data
                    map_roicom = mapd[
                        round(roicom[0]), round(roicom[1]), round(roicom[2])
                    ]

                # *** Continue calculating the overlap of the ROI and the map
                # Only if COM was not enabled, or when it was and the COM
                # overlaps with the map data
                if not self.args.com or map_roicom == 1:

                    # Use the ROI as a mask for the map.
                    map_masked = np.multiply(mapd, roibin)

                    # Count the number of voxels in this masked image
                    map_masked_size = np.sum(map_masked)

                    # Calculate what % this is of the original ROI size
                    map_masked_perc = (map_masked_size / roisize) * 100

                    # If this % passes the user threshold select this ROI
                    if map_masked_perc >= thr:
                        select_roi = 1

                    else:
                        select_roi = 0

                # ...in all other cases, don't select this ROI
                else:
                    select_roi = 0
                    map_masked_size = 0
                    map_masked_perc = 0

                # Return output
                return([select_roi, roisize, map_masked_size, map_masked_perc])
            # Add function to main the object
            self.test_roi = test_roi

            # ** Output array to store which ROIs were selected
            self.output1 = np.empty((np.size(self.atlas_data.roinums), 5))

            # ** Function to run the overlap-test function for all ROIs
            def run_test_roi(roinums, atld, mapd, thr):

                # *** Start row index
                outrow = 0

                # *** Loop over ROIs
                for i in roinums:

                    # **** Store the ROI number in the output matrix
                    self.output1[outrow, 0] = i

                    # **** Test for overlap
                    select_roi = self.test_roi(i, atld, mapd, thr)

                    # **** Store output in the output matrix
                    self.output1[outrow, 1:5] = select_roi

                    # **** Continue to the next row on the output matrix
                    outrow = outrow + 1
            # Add function to the main object
            self.run_test_roi = run_test_roi

            # ** Run the funtion to do all overlap testing
            self.run_test_roi(
                self.atlas_data.roinums,
                self.atlas_data.imgrdata.get_fdata(),
                self.map_data.binimg,
                self.args.othr1
            )

            # ** Save output matrix to output file in output folder
            ofile = self.args.outdir + '/ROI_selectionR.csv'
            np.savetxt(
                ofile, self.output1, delimiter=',', fmt='%i,%i,%i,%i,%3.2f'
            )

            # ** Add header
            # I don't want to use Pandas to include a header, so I will just add
            # a header to the text file.
            with open(ofile) as fp:
                data = fp.read()

            header = "ROI value,Select R,ROI size (v),Overlap (v),Overlap (%)"

            with open(ofile, 'w') as fp:
                fp.write(header + "\n" + data)

        # ** -------------------------------------------------------------------
        # ** PART 2: Select ROIs that contain a sufficiently large portion of
        # ** the input image
        # ** -------------------------------------------------------------------
        if args.parts == "both":

            # Select all ROIs that individually comprise a certain
            # (user-defined) part of a free-floating cluster of the input map
            # image. This can be helpful if the map contains small clusters that
            # are of interest, but which are not big enough to fill an entire
            # ROI for the user-defined threshold percentage.

            # ** Announce
            print("PART 2: Selecting ROIs that contain a sufficiently large "
                  "portion of the input image...")

            # ** Convert the input map into individual clusters
            # A cluster is defined as a cluster of ones surrounded by zeros in
            # the user input map, after binarization based on the user set
            # threshold.
            class map_cluster:
                def __init__(self, maptr, outdir):

                    # *** Output file name
                    self.outfile = outdir + '/map_tc.nii.gz'

                    # *** Extract clustes
                    self.cl = fsl.Cluster()
                    self.cl.inputs.threshold = 1
                    self.cl.inputs.in_file = maptr
                    self.cl.inputs.out_index_file = self.outfile
                    self.cl.terminal_output = 'none'
                    self.results = self.cl.run()

                    # *** List the number of clusters
                    self.imgc = nb.load(self.outfile)
                    self.cnums = np.unique(self.imgc.get_fdata())
                    # Remove zero, because that is not a label
                    self.cnums = self.cnums[self.cnums != 0]
            # Call: Convert the input map to individual clusters
            self.map_cluster = map_cluster(
                self.map_data.outfile, self.args.outdir
            )

            # ** Function to test if a cluster fits inside an ROI above threshld
            def test_cluster(mapc, c, atld, i, thr):

                # *** Extract cluster
                clubin = np.where(mapc == c, 1, 0)

                # *** Size of the cluster in voxels
                clusize = np.sum(clubin)
                # Number of voxels at minimal overlap
                clusize_t = (clusize / 100) * thr

                # *** Test if it is even possible for the cluster to fit the ROI
                # If the number of voxels of the cluster at the percentage
                # threshold is bigger than the ROI, it will never fit for that
                # minimum percentage inside the ROI.
                # E.g., if the cluster=10 voxels, and the ROI=2 voxels, the
                # cluster (even at a trheshold of 50%) would not fit inside the
                # ROI. We can quickly test this, because we have the size of the
                # cluster and the size of the ROI in our output matrix.

                # *** Look up ROI size
                roisize = self.output1[self.output1[:, 0] == i][0, 2]

                # *** Only continue if the ROI size is bigger than the cluster
                if roisize > clusize_t:

                    # *** Extract ROI (binary)
                    roibin = np.where(atld == i, 1, 0)

                    # Use the ROI as a mask for the cluster map.
                    map_masked = np.multiply(clubin, roibin)

                    # Count the number of voxels in this masked image
                    map_masked_size = np.sum(map_masked)

                    # Calculate what % this is of the original cluster size
                    map_masked_perc = (map_masked_size / clusize) * 100

                    # If this % passes the user threshold select this ROI
                    if map_masked_perc >= thr:
                        select_roi = 1

                    else:
                        select_roi = 0

                # ...in all other cases, don't select this ROI
                else:
                    select_roi = 0
                    map_masked_size = 0
                    map_masked_perc = 0

                # Return output
                return([select_roi, roisize, clusize, clusize_t,
                        map_masked_size, map_masked_perc])
            # Add function to main the object
            self.test_cluster = test_cluster

            # ** Output array to store which ROIs were selected
            self.output2 = np.empty((
                np.size(
                    self.atlas_data.roinums
                ) * np.size(
                    self.map_cluster.cnums
                ),
                8
            ))

            # ** Function to run the overlap function for all clusters and ROIs
            def run_test_clu(clunums, roinums, mapc, atld, thr):

                # *** Start row index
                outrow = 0

                # *** Loop over clusters
                for c in clunums:

                    # **** Loop over ROIs
                    for i in roinums:

                        # ***** Store the Cluster and ROI number in the matrix
                        self.output2[outrow, 0] = c
                        self.output2[outrow, 1] = i

                        # ***** Test for overlap
                        select_roi = self.test_cluster(mapc, c, atld, i, thr)

                        # ***** Store output in the output matrix
                        self.output2[outrow, 2:9] = select_roi

                        # ***** Continue to the next row on the output matrix
                        outrow = outrow + 1
            # Add function to the main object
            self.run_test_clu = run_test_clu

            # ** Run the funtion to do all overlap testing
            self.run_test_clu(
                self.map_cluster.cnums,
                self.atlas_data.roinums,
                self.map_cluster.imgc.get_fdata(),
                self.atlas_data.imgrdata.get_fdata(),
                self.args.othr2
            )

            # ** Save output matrix to output file in output folder
            ofile = self.args.outdir + '/ROI_selectionC.csv'
            np.savetxt(
                ofile, self.output2, delimiter=',',
                fmt='%i,%i,%i,%i,%i,%i,%i,%3.2f'
            )

            # ** Add header
            # I don't want to use Pandas to include a header, so I will just add
            # a header to the text file.
            with open(ofile) as fp:
                data = fp.read()

            header = 'Cluster,ROI value,Select R,ROI size (v),' \
                'Cluster size (v), Cluster size thr (v),Overlap (v),' \
                'Overlap (%)'

            with open(ofile, 'w') as fp:
                fp.write(header + "\n" + data)

        # ** -------------------------------------------------------------------
        # ** PART 3: Create an output image with all clusters
        # ** -------------------------------------------------------------------

        # ** Announce
        print("Creating an output nifti image with all selected clusters...")
        print("\nSelected ROIs:")

        # ** Function to create output image
        # That is, all selected ROIs with their original values
        def build_outimg(omatrix, atld):

            # ** Create empty output image
            outimg = np.zeros(self.map_data.binimg.shape)

            # ** Grab list of seleced ROIs
            selection1 = np.array([])
            selection2 = np.array([])

            # *** Selected ROIs for Part 1
            if args.parts == "p1" or args.parts == "both":
                selection1 = self.output1[self.output1[:, 1] == 1]
                selection1 = selection1[:, 0]
                print("Part 1: [n=" + str(len(selection1)) + "]:  " + " ".join(
                    map(str, selection1)).replace(".0", " "))

            # *** Selected ROIs for Part 2
            if args.parts == "both":
                selection2 = self.output2[self.output2[:, 2] == 1]
                selection2 = selection2[:, 1]
                print("Part 2: [n=" + str(len(selection2)) + "]:  " + " ".join(
                    map(str, selection2)).replace(".0", " "))

            # *** Combine selections
            selection = np.concatenate((selection1, selection2))
            selection = np.unique(selection)
            if args.parts == "both":
                print("Total:  [n=" + str(len(selection)) + "]:  " + " ".join(
                    map(str, selection)).replace(".0", " "))

            # ** Loop over selected ROIs and add them together
            for ROI in selection:

                # *** Extract ROI (original value)
                roiimg = np.where(atld == ROI, atld, 0)

                # *** Add up
                outimg = outimg + roiimg

            return(outimg)
        # Add function to the main object
        self.build_outimg = build_outimg

        # ** Save output image to output file in output folder
        ofile = self.args.outdir + '/ROI_selection.nii.gz'
        odata = nb.Nifti1Image(
            self.build_outimg(
                self.output1,
                self.atlas_data.imgrdata.get_fdata()
            ),
            self.map_data.img.affine,
            self.map_data.img.header
        )
        nb.save(odata, ofile)


# * Run code
m2a = map2atlas()
