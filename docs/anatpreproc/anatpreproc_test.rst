.. AUTO-GENERATED FILE -- DO NOT EDIT!

.. _example_anatpreproc_test:


Anatomical Workflow Sanity Checks
================================

**Class anatPreprocTest**
-------------------------


::
  
  class anatPreprocTest(object):
  

This class defines all the quantitative checks done to the inputs and outputs of
Anatomical Preprocessing

::
  
  def __init__(self, preproc, base_dir, input_anat):

Constructor to initialize the workflow

**Parameters**

 *preproc* : (Object)
              Anatomical workflow object
 *base_dir* : (String)
             path of output workflow directory
 *input_anat* : (file)
               Input image to the anatomical workflow


  

Setup workflow object
---------------------

::
  
  def setup(self):

Set up the Workflow and run
This method is run once before _each_ test method is executed

**Example**

>>> self.preproc.inputs.anat = self.input_anat
>>> self.preproc.base_dir = self.base_dir
>>> self.preproc.run()


  

Delete workflow  object
-----------------------

::
  
  def teardown(self):

Delete The Workflow Object
This method is run once after _each_ test method is executed

  

Test to verify inputs
---------------------

::
  
  def inputs_test(self):
  assert False

Method to check if the input file
is T1 image and is in right format

**Returns**

  *TRUE* : (boolean)
           if all the tests pass
  *Error* : (Exception)
            raise if any of the tests fail

**Tests**

  - input_anat image should be a nifti file with extension nii.gz or nii
  - input_anat image should be a 3D image, i.e only one volume

  

Test to verify deoblique Image
------------------------------

::
  
  def anat_deoblique_test(self, deoblique_img):
  assert False

method to check if the input file
is deobliqued correctly or not

**Parameters**

*deoblique_img* : (nifti file)
                  De-obliqued mprage file

**Returns**

  Same as above method

**Tests**

  - Compare the headers of input_anat and deoblique_img,
    the voxel offset values will be unequal, if the original
    image is obliqued. All other values should remain unchanged.
  - Compare the Image data of input_anat and deoblique_img, they
    should remain same


  

Test to verify reorientation
----------------------------

::
  
  def anat_reorient_test(self, standard_img, rpi_img):
  assert False

method to check if the output reorient is coorectly
in RPI orientation or not

**Parameters**

 *standard_img* : (nifti file)
                  Standard T1 RPI image in MNI space
 *rpi_img* : (nifti file)
             RPI output mprage file

**Returns**

  Same as above method

**Tests**

 - Compute Spatial Correlation between standard_img and rpi_img
   For this, first convert the RPI output image into MNI space
   and if necessary resample it to match the voxel dimensions
   of the standard image
 - Compare the Image data of rpi_img and input_anat. It should be
   same.


  

Test to verify skulStrip image with normalized/scaled intensities
-----------------------------------------------------------------

::
  
  def anat_skullstrip_test(self, skullstrip_img):
  assert False  

method to check if the skull stripped image is
correct or not

**Parameters**

 *skullstrip_img* : (nifti file)
                    Skullstrip output image with normalized/scaled intensities

**Returns**

  Same as above method

**Tests**

  - Since afni scales/normalizes the intensity values
    its hard to test it.Can be tested in the next step

  

Test to verify skullstrip image with original intensity values
--------------------------------------------------------------

::
  
  def anat_brain_test(self, rpi_img, brian_img):
  assert False

method to check if the input file
is skull stripped and the intensity values are unchanged

**Parameters**

 *rpi_img* : (nifti file)
             RPI output mprage file
 *brian_img* : (nifti file)
               Skull stripped Brain only file

**Returns**

  Same as above method

**Tests**

 - Subtract (infile_a - infile_b), this should return a matrix
   with all zero values in brain area and non zero values around
   the edge. From the origin, choose a sphere of resonable diameter
   and check the intensity values should be zero. Then check the
   edges for non-zero values.

