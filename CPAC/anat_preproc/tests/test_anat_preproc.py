import os
import sys
from  .. import anat_preproc
from nose.tools import *
import nibabel as nib
import numpy as np

class TestAnatPreproc:


    def __init__(self):

        """
            Initialize and run the the anat_preproc workflow

            Populate the node-name : node_output dictionary using the workflow object.
            This dictionary serves the outputs of each of the nodes in the workflow to all the tests that need them.

            Parameters
            ----------

            self

            Returns
            -------

            None


        """


        self.preproc = anat_preproc.create_anat_preproc()
        self.input_anat = os.path.abspath('/usr/share/data/fsl-mni152-templates/MNI152_T1_2mm.nii.gz')
        self.input_anat_known_brain = os.path.abspath('/usr/share/data/fsl-mni152-templates/MNI152_T1_2mm_brain.nii.gz')
        self.base_dir = os.path.abspath('/home/data/Projects/Pipelines_testing/Dickstein_CPAC_working_dir/')
        self.preproc.inputs.inputspec.anat = self.input_anat
        self.preproc.base_dir = self.base_dir
        self.anat_execution_graph = self.preproc.run()

        self.node_names_map = {}

        it = self.anat_execution_graph.nodes_iter()

        for node in it:

            name = node.fullname

            node_output = node.result.outputs.out_file

            name = name.split('.')[1]

            #print node.result.outputs.out_file, ' ', node.fullname

            if not (name in self.node_names_map):
                self.node_names_map[name] = node_output


#    def setUp(self):



    def tearDown(self):
        """
            Deletes the workflow object
            
        """

        del self.preproc

    def test_inputs(self):

        """

            Checks the input anatomical image.

            Parameters
            ----------

            self

            Returns
            -------

            None

            Notes
            -----

            The tests are 

            - Check if the input image ends with .nii or .nii.gz
            - Check if the input image is stored in 3D only 


        """


        anat_img = self.input_anat

        if not (anat_img.endswith('.nii.gz') or anat_img.endswith('.nii')):

            raise 'Input image name does not end with .nii.gz or .nii'


        try:

            img = nib.load(anat_img)
            dims = img.get_shape()

            assert len(dims) == 3

        except:

            raise 'input T1 image is not 3D'



    def test_anat_deoblique(self):

        """
            Check the output from anat_deoblique node

            Parameters
            ----------

            self

            Returns
            -------

            None


            Notes
            -----

            3drefit command is supposed to overwrite the transformation matrix with the cardinal matrix in the header
            In addition to this it also modifies a lot of other fields(aside from voxel_offset it sets other header fields
            to zero) . The number of fields varies depending on how the nifti file was obtained from the source dicom files.


            So, The check tests the image Data of the output of anat_deoblique node against the input anatomical image.
            Both should be the same

        """


        if not ('anat_deoblique' in self.node_names_map):
            assert False

        else:
            deobliqued_anat = self.node_names_map['anat_deoblique']

            de_img_data = nib.load(deobliqued_anat).get_data()
            orig_image_data = nib.load(self.input_anat).get_data()


            de_img_data = de_img_data.flatten()
            orig_image_data = orig_image_data.flatten()

            np.testing.assert_equal(de_img_data, orig_image_data)

#            de_image_header = nib.load(deobliqued_anat).get_header()
#
#            de_vox_offset = de_image_header.get_data_offset()
#
#            orig_image_header = nib.load(self.input_anat).get_header()
#
#            orig_vox_offset = orig_image_header.get_data_offset()
#
#
#            if orig_vox_offset == de_vox_offset:
#
#                print 'Its not possible for deobliqued T1 and Original T1 voxels offsets are same'
#                raise
#
#            else:
#                    not_match_count = 0
#                    header_values_de = de_image_header.values()
#                    header_values_input = orig_image_header.values()
#
#                    for i in range(0, len(header_values_input)):
#
#                        a = (header_values_de[i] == header_values_input[i])
#                        if not (a.all() == True):
#
#                            not_match_count += 1
#                    print 'not_match_count: ', not_match_count
#                    #assert not_match_count == 1



    def test_anat_reorient(self):

        """
            Checks the output of anat_reorient node.

            Parameters
            ----------

            self

            Returns
            -------

            None


            Notes
            -----

            Checks the s-form matrix in the header of the output file.
            For the output image to be in RPI , The first diagonal element must be negative
            The rest of the diagonal elements must be positive. The non diagonal elements must be
            zero, the offset elements can be anything.

        """



        if not ('anat_reorient' in self.node_names_map):
            assert False

        else:
            anat_reorient = self.node_names_map['anat_reorient']
            anat_reorient_sform = nib.load(anat_reorient).get_sform()

            #print 'sform: ', anat_reorient_sform

            for i in range(0, 4):

                for j in range(0, 4):

                    if i == j:
                        if i == 0:

                            assert int(anat_reorient_sform[i][i]) < 0
                        else:
                            assert int(anat_reorient_sform[i][i]) > 0

                    else:

                        if not (j == 3):

                            assert int(anat_reorient_sform[i][j]) == 0



    def test_anat_skullstrip(self):

        """

            Tests the output of anat_skullstrip node in the workflow

            Parameters
            ----------

            self

            Returns
            -------

            None

            Notes
            -----

            Compares the output image of anat_skullstrip node with the known skullstripped image.
            Since the intensity values do not correspond to the input 

        """


        if not ('anat_skullstrip' in self.node_names_map):

            assert False

        else:

            def binarize(data):

                data[np.flatnonzero(data)] = 1
                return data

            anat_skullstrip = self.node_names_map['anat_skullstrip']

            skullstrip_data = nib.load(anat_skullstrip).get_data()
            known_brain_data = nib.load(self.input_anat_known_brain).get_data()

            bin_skullstrip = binarize(skullstrip_data.flatten())
            bin_brain = binarize(known_brain_data.flatten())

            correlation = np.corrcoef(bin_skullstrip, bin_brain)

            #print 'correlation: ', correlation
            assert correlation[0, 1] >= 0.95



    def test_anat_brain(self):

        """
            Compares the SkullStripped Anatomical Image(after its passed through step function) with the known skullstripped image available

            Parameters
            ----------

            self

            Returns
            -------

            Nothing
        """

        if not ('anat_brain_only' in self.node_names_map):

            assert False

        else:


            anat_brain = self.node_names_map['anat_brain_only']

            brain_data = nib.load(anat_brain).get_data()
            known_brain_data = nib.load(self.input_anat_known_brain).get_data()


            correlation = np.corrcoef(brain_data.flatten(), known_brain_data.flatten())

            #print 'correlation: ', correlation

            assert correlation[0, 1] >= 0.97
