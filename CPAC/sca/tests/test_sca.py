import sys
#import e_afni
import os
import subprocess
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
#from utils import *

class TestSCA(object):

    def __init__(self, rest_res_filt, ref, standard, fwhm, seed_list, rest_mask2standard, premat, postmat, extraction_space):

        """
        
        Initialize Inputs
        Constructor call to initialize the inputs to sca_preproc workflow
        
        Parameters
        ----------
         
        seed_list : a list of existing nifti files
               A list of seeds/ ROI iin MNI space.
        
        extraction_space : (a string)
            Options are 'mni' or 'native'.
               Extract time series from ROI in MNI space or extract it in subjects native space.
        
        rest_res_filt : (an existing nifti file) 
            Band passed Image with Global Signal , white matter, csf and motion regression. Recommended bandpass filter (0.001,0.1) )
               (rest_res_bandpassed.nii.gz)
        
        ref : (an existing nifti file)
            When Registering from MNI space to native space use the mean functional image in native space example_func.nii.gz is used.
        
        standard : (an existing nifti file)
            When registering from native to MNI MNI152_T1_STANDARD_RES.nii.gz is used(target space).
        
        fwhm : (A list of floating point numbers)
            For spatial smoothing the Z-transformed correlations in MNI space.
               Generally the value of this parameter is 1.5 or 2 times the voxel size of the input Image.
        
        Returns
        -------
        
        Nothing            
        

        """

        self.rest_res_filt = rest_res_filt
        self.ref = ref
        self.standard = standard
        self.fwhm = fwhm
        self.premat = premat
        self.postmat = postmat
        self.seed_list = seed_list
        self.rest_mask2standard = rest_mask2standard
        self.extraction_space = extraction_space
        self.sca_preproc = self.setup_workflow()


    def setup_workflow(self):
        """
        Set up the workflow object by initializing it

        Parameters
        ----------

        self

        Returns
        -------

        None

        """

        self.sca_preproc = sca_preproc(self.extraction_space)


    def teardown_workflow(self):
        """

        Tear down the workflow by deleting the workflow object

        Parameters
        ----------

        self

        Returns
        -------

        None

        """
        del self.sca_preproc

    def test_inputs(self):
        """

        Test the workflow inputs
 
        Parameters
        ----------

        self

        Returns
        -------
 
        Generates exceptions if any of the input tests fails 

        Notes
        -----
    
        seed_list :

             Verify that inputs seeds lie inside the MNI brain(Not clear on how to verify it without the user supplied co-ordinates.
            The ROI can be in native space but indicate MNI space in the header).
            The input seed files are nifti files(.nii or .nii.gz).
             All the seeds have the same resolution.
            
        extraction_space :

            Must be a string uppercase or lower case
            Options are 'mni' or 'native'
            
        rest_res_filt :

            Test if the file is nifti file. Since band pass filtering can be turn on and off, we should decide if we need to test for temporal filtering
            Test if the file is nuisance signal regressed according to the nuisance signal corrections specified(This test would duplicate testing efforts
            between the nuisance signal regression workflow tests)
            
        ref :

            Test if the file is nifti file.
            Test the co-ordinate space of the image; Compute spatial correlation with the Standard MNI image in the given resolution.
            If the spatial correlation is 1 then the reference file is verified to be Standard MNI image; return 1.
            Otherwise compute correlation between ref image and example_func nifti image generated from functional preprocessing workflow.
            If the spatial correlation is 1 or epsilon close to 1 then return 1 else return 0.
        
        fwhm :
        
            Test if this is a list of floating point numbers.
        
        
        """

        assert False

    def test_output(self):
        """
        
        Run the sca_preproc workflow and test all the outputs 
       
        Parameters 
        ----------

        self    

        Returns
        -------
            
        generate exceptions for each outputs that fail the test

        Notes
        -----

        correlations_verification : (a boolean value)
            Reports if correlations command works as expected.
                
        z_trans_correlations_verication : (a boolean value)
            Reports if Image normalization works as expected.
        
        z_2standard_verification : (a boolean value)
            Reports if registration to template space works as expected.
                
        z_2standard_FWHM_verification : (a boolean value)
            Reports if spatial smoothing on z_2standard image works as expected
        
        Detailed_Description
        
        - Verify Correlations :
        
            On Real Data -
        
            Extract mean TimeSeries using the ROI from rest_res_filt nifti image
            Compute the map of Voxel wise correlations using nibabel and numpy
            Compare the correlations calculated manually against the ones in the input Correlations image(correlations nifti image)
            Since Comparison is between floating points keep a + -  epsilon slack in comparison; epsilon ~= 0.005
            Verify all values are between 0 & 1
            If both these test evaluate to True then return True else Return False
        
            On Simulation Data -
        
            Not sure

        - Verify Image Normalization :
        
            On Real Data -

            Use the voxel wise Correlation map from Correlations test to manually compute fisher-Z normalization scores
            Compare the calculated correlation map values against the values from Z normalized input image
            Check if the variance in the z scores is approx constant and that the correlation among values calculated in the prev two steps is epsilon closer to one 
        
            On Simulation Data -
        
            Not sure
        
        - Verify Image Registration :
        
            On Real Data -
            
            Compute the spatial correlation of the registered Image(z_2standard) with the reference nifti image (ref)
            Compute percentage overlap between registered image and the reference image and
            between the input image(z_trans_correlations) and the reference image
            The correlation values obtained when registered image is used should be higher than the when the input image is used (Need to test).
            If the spatial correlation and percentage overlaps are above their respective thresholds then return True else False
            (TODO: Decide before hand on the values of thresholds, a way to do that could be an average of correlation values of subjects on a big dataset)
        
        
            On Simulation Data -
            
            Not sure
            
            
        - Verify Guassian Kernel Smoothing :
        
            On Real Data -
            
            Verify that the intensiy distribution in smoothed image is still centered at zero and have unit variance
            Load the z_2standard input image in python. Form the (2-D ) Guassian kernel matrix( TODO: determine the size of the gaussian matrix)
            Perform convolution of the loaded input image and the gaussian matrix to get the guassian smoothed image
            Perform mean filtering (substitute the intensity of the voxel with the average intesity of the neighbouring voxels)
            Compare the image yeilded from the previous step with the smoothed image input image(z_2standard_FWHM)
            if the values are +- epsilon away then return True else return False
        
            On Simulation Data -
        
            Not sure
        
        
        """

        assert False



    def test_warp_to_native(self):
        """
        Checks if warping to native space works as expected
        Set up the workflow object by initializing it

        Parameters
        ----------

        self

        Returns
        -------

        generates an exception if warp_to_native_test fails

        Notes
        -----

        On Real Data -
        
        Compute the spatial correlation of the registered Image(outputs of the warp_to_native node) with the reference nifti image (ref)
        Compute percentage overlap between registered image and the reference image and
        between the input image(Seeds in MNI space) and the reference image
        The correlation values obtained when registered image is used should be higher than the when the input image is used (Need to test).
        
        On Synthetic Data -
        
        - Not sure


        """

        assert False


    def test_time_series(self):
        """
        checks if extraction of timeseries works as expected
        
        Parameters
        ----------
        
        self
        
        Returns
        -------
        
        generates an exception if time_series_test fails

        Notes
        -----

        On Real Data -
        
        CASE_extraction_space_MNI: Extract the mean time series from the rest_res_filt in MNI space(do it for all the seeds in MNI space)
        and find the correlation with the output of the timeseries node.
        If the correlation is epsilon close to 1 then the data passes the test.
        
        CASE_extraction_space_Native: Extract the mean time series from the rest_res_filt (warp all the seeds to Native space first)
        and find the correlation with the output of the timeseries node.
        If the correlation is epsilon close to 1 then the data passes the test.
        
        On Simulation Data -
        
        - Not sure 
        """

        assert False


    def test_print_timeseries_to_file(self):
        """
        checks if outputing timeseries list to a file works as expected
        
        Parameters
        ----------
        
        self
        
        Returns
        -------
        
        generates an exception if  print_timeseries_to_file_test fails
        
        Notes
        -----
        
        On Real Data -
        
        Read the timeseries from the timeseries file (output of print_timeseries_to_file node).
        Correlate with output of timeseries node , The test is verified if the correlation is epsilon close to 1
        
        On Simulation Data -
        
        Not sure 
        
        """

        assert False


    def test_warp_filt(self):


        """
        checks if warping from functional space to MNI space works as expected
        
        Parameters
        ----------
        
        self
        
        Returns
        -------
        
        generates an exception if warp_filt_test fails
        
        Notes
        -----
        
        On Real Data -
        
        Compute the spatial correlation of the registered Image(output of warp_filt node) with the reference nifti image (standard)
        Compute percentage overlap between registered image and the reference image and
        between the input image(rest_res_filt) and the reference image
        The correlation values obtained when registered image is used should be higher than the when the input image is used (Need to test).
        Also check the level of noise introduced in the bandpassed frequency band. If noise significant report bad registration( Need to define what significant is)
        
        On Simulation Data -
        
        Not Sure
        
        """

        assert False


    def test_z_trans(self):

        """
        checks if fisher Z transformation works as expected
        
        Parameters
        ----------
        
        self
        
        Returns
        -------
        
        generates an exception if z_trans_test fails
        
        Notes
        -----
        
        On Real Data -
        
        Measure the variance in z scores in a subject and between subjects as the correlation values vary, is approximately constant
        
        On Simulation Data -
        
        Not sure    
        
        """
        assert False


    def test_warp_to_standard(self):

        """
        checks if warping from functional space to MNI space works as expected
        
        
        Parameters
        ----------
        
        self
        
        Returns
        -------
        
        generates an exception if warp_to_standard_test fails
        
        Notes
        -----
        
        On Real Data -
        
        Compute the spatial correlation of the registered Image(output of warp_to_standard node) with the reference nifti image (standard)
        Compute percentage overlap between registered image and the reference image and
        between the (z_trans_correlations) and the reference image
        The correlation values obtained when registered image is used should be higher than the when the input image is used (Need to test).
        Also verify that variation in intesity values is approximately constant(does not shoot up as the intensity approaches 1 and vice versa as it approaches -1)
        in a subject and between subjects 
        
        On Simulation Data -
        
        Not Sure
        
        """
        assert False


    def test_corr(self):

        """
        tests if computing correlations works as expected
        
        Parameters
        ----------
        
        self
        
        Returns
        -------
        
        generates an exception if corr_test fails
        
        Notes
        -----
        
        On Real Data -
        
        Extract mean TimeSeries using the ROI from rest_res_filt nifti image(use the seed and nifi image in appropriate co-ordinate space for extraction_space = 'mni' & 'native')
        Compute the map of Voxel wise correlations using nibabel and numpy
        Compare the correlations calculated manually against the ones in the input Correlations image(correlations nifti image)
        Since Comparison is between floating points keep a + -  epsilon slack in comparison; epsilon ~= 0.005
        Verify all values are between 0 & 1
        If both these test evaluate to True then return True else Return False
        
        On Simulation Data -
        
        Not sure
        
        
        """
        assert False


    def test_smooth_mni(self):

        """
        test if spatial smoothing works as expected
        
        Parameters
        ----------
        
        self
        
        Returns
        -------
        
        generates an exception if smooth_mni_test fails
        
        Notes
        -----
        
        On Real Data -
        
        Verify that the intensiy distribution in smoothed image is still centered at zero and have unit variance
        Load the z_2standard input image in python. Form the (2-D ) Guassian kernel matrix( TODO: determine the size of the gaussian matrix)
        Perform convolution of the loaded input image and the gaussian matrix to get the guassian smoothed image
        Perform mean filtering (substitute the intensity of the voxel with the average intesity of the neighbouring voxels)
        Compare the image yeilded from the previous step with the smoothed image input image(z_2standard_FWHM)
        if the values are +- epsilon away then return True else return False
        
        On Simulation Data -
        
        Not sure
        
        
        """
        assert False
