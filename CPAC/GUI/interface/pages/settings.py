import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import os
import pkg_resources as p


class Settings(wx.html.HtmlWindow):

    def __init__(self, parent, counter=0):
        wx.html.HtmlWindow.__init__(
            self, parent, style=wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()

        self.counter = counter

        self.LoadFile(p.resource_filename('CPAC',
                                          'GUI/resources/html/compute_config.html'))

    def get_counter(self):
        return self.counter


class ComputerSettings(wx.ScrolledWindow):

    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)
        self.counter = counter
        
        fsl = ""
        if os.environ.get('FSLDIR'):
            fsl = os.environ['FSLDIR']

        self.page = GenericClass(self, "Computer Settings")
        self.page.add(label="Run CPAC on a Cluster/Grid ",
                      control=control.CHOICE_BOX,
                      name='runOnGrid',
                      type=dtype.BOOL,
                      comment="Select False if you intend to run CPAC on a "
                              "single machine.\n\nIf set to True, CPAC will "
                              "attempt to submit jobs through the job "
                              "scheduler / resource manager selected below.",
                      values=["On", "Off"],
                      wkf_switch=True)

        self.page.add(label="FSL Path ",
                      control=control.DIR_COMBO_BOX,
                      name='FSLDIR',
                      type=dtype.STR,
                      values=fsl,
                      comment="Full path to the FSL version to be used by "
                              "CPAC.\n\nIf you have specified an FSL path "
                              "in your .bashrc file, this path will be set "
                              "automatically.")

        self.page.add(label="Job Scheduler / Resource Manager ",
                      control=control.CHOICE_BOX,
                      name='resourceManager',
                      type=dtype.STR,
                      values=["SGE", "PBS", "SLURM"],
                      comment="Sun Grid Engine (SGE), "\
                              "Portable Batch System (PBS), or "
                              "Simple Linux Utility for Resource Management (SLURM)."\
                              "\n\nOnly applies if you are running on a grid or compute cluster.")

        self.page.add(label='SGE Parallel Environment ',
                      control=control.TEXT_BOX,
                      name='parallelEnvironment',
                      type=dtype.STR,
                      comment='SGE Parallel Environment to use when running CPAC.\n\nOnly applies when you are running on a grid or compute cluster using SGE.',
                      values='cpac')

        self.page.add(label='SGE Queue ',
                      control=control.TEXT_BOX,
                      name='queue',
                      type=dtype.STR,
                      comment='SGE Queue to use when running CPAC.\n\nOnly applies when you are running on a grid or compute cluster using SGE.',
                      values='all.q')

        self.page.add(label="Maximum Memory Per Participant (GB) ",
                      control=control.INT_CTRL,
                      name='maximumMemoryPerParticipant',
                      type=dtype.NUM,
                      comment="The maximum amount of memory each " \
                              "participant's workflow can allocate. Use " \
                              "this to place an upper bound of memory " \
                              "usage. Warning: 'Memory Per Participant' " \
                              "multiplied by 'Number of Participants to " \
                              "Run Simultaneously' must not be more than " \
                              "the total amount of RAM. Conversely, using " \
                              "too little RAM can impede the speed of a " \
                              "pipeline run. It is recommended that you set "\
                              "this to a value that when multiplied by " \
                              "'Number of Participants to Run " \
                              "Simultaneously' is as much RAM you can " \
                              "safely allocate.",
                      values=1)

        self.page.add(label="Maximum Number of Cores Per Participant ",
                      control=control.INT_CTRL,
                      name='maxCoresPerParticipant',
                      type=dtype.NUM,
                      comment="The maximum amount of cores (on a single " \
                              "machine) or slots on a node (on a " \
                              "cluster/grid) to allocate per participant. " \
                              "Setting this above 1 will parallelize each " \
                              "participant's workflow where possible. If " \
                              "you wish to dedicate multiple cores to " \
                              "ANTS-based anatomical registration (below), " \
                              "this value must be equal or higher than the " \
                              "amount of cores provided to ANTS. The " \
                              "maximum number of cores your run can " \
                              "possibly employ will be this setting " \
                              "multiplied by the number of participants " \
                              "set to run in parallel (the 'Number of" \
                              "Participants to Run Simultaneously' " \
                              "setting).",
                      values=1)

        self.page.add(label="Number of Participants to Run Simultaneously ",
                      control=control.INT_CTRL,
                      name='numParticipantsAtOnce',
                      type=dtype.NUM,
                      comment="The number of participant workflows to run " \
                              "at the same time. The maximum number of " \
                              "cores your run can possibly employ will be " \
                              "this setting multiplied by the number of " \
                              "cores dedicated to each participant (the " \
                              "'Maximum Number of Cores Per Participant' " \
                              "setting).",
                      values=1)

        self.page.add(label="Number of Cores for Anatomical Registration " \
                            "(ANTS only) ",
                      control=control.INT_CTRL,
                      name='num_ants_threads',
                      type=dtype.NUM,
                      comment="The number of cores to allocate to " \
                              "ANTS-based anatomical registration per " \
                              "participant. Multiple cores can greatly " \
                              "speed up this preprocessing step. This " \
                              "number cannot be greater than the number of " \
                              "cores per participant.",
                      values=1)

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter


class DirectorySettings(wx.ScrolledWindow):

    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)

        self.page = GenericClass(self, "Output Settings")
        self.counter = counter
        
        self.page.add(label="Pipeline Name ",
                      control=control.DIR_COMBO_BOX,
                      name='pipelineName',
                      type=dtype.STR,
                      comment="Name for this pipeline configuration - useful for identification.",
                      validation_req=False)

        self.page.add(label="Working Directory ",
                      control=control.DIR_COMBO_BOX,
                      name='workingDirectory',
                      type=dtype.STR,
                      comment="Directory where CPAC should store temporary and intermediate files.",
                      validation_req=False)

        self.page.add(label="Crash Log Directory ",
                      control=control.DIR_COMBO_BOX,
                      name='crashLogDirectory',
                      type=dtype.STR,
                      comment="Directory where CPAC should write crash logs.",
                      validation_req=False)

        self.page.add(label="Log Directory ",
                      control=control.DIR_COMBO_BOX,
                      name='logDirectory',
                      type=dtype.STR,
                      comment="Directory where CPAC should place run logs.",
                      validation_req=False)

        self.page.add(label="Output Directory ",
                      control=control.DIR_COMBO_BOX,
                      name='outputDirectory',
                      type=dtype.STR,
                      comment="Directory where CPAC should place processed "
                              "data.",
                      validation_req=False)

        self.page.add(label="AWS Output Bucket Credentials (optional) ",
                      control=control.COMBO_BOX,
                      name='awsOutputBucketCredentials',
                      type=dtype.STR,
                      comment="If setting the \'Output Directory\' to an S3 "
                              "bucket, insert the path to your AWS "
                              "credentials file here.",
                      validation_req=False)

        self.page.add(label="S3 Encryption ",
                      control=control.CHOICE_BOX,
                      name='s3Encryption',
                      type=dtype.LSTR,
                      comment="Enable server-side 256-AES encryption on data "
                              "to the S3 bucket",
                      values=["On", "Off"])

        self.page.add(label="Write Extra Functional Outputs ",
                      control=control.CHOICE_BOX,
                      name='write_func_outputs',
                      type=dtype.LSTR,
                      comment="Include extra versions and intermediate steps "
                              "of functional preprocessing in the output "
                              "directory.",
                      values=["Off", "On"])

        self.page.add(label="Write Debugging Outputs ",
                      control=control.CHOICE_BOX,
                      name='write_debugging_outputs',
                      type=dtype.LSTR,
                      comment="Include extra outputs in the output "
                              "directory that may be of interest when more "
                              "information is needed.",
                      values=["Off", "On"])

        self.page.add(label="Enable Quality Control Interface ",
                      control=control.CHOICE_BOX,
                      name='generateQualityControlImages',
                      type=dtype.LSTR,
                      comment="Generate quality control pages containing "
                              "preprocessing and derivative outputs.",
                      values=["On", "Off"])

        self.page.add(label="Remove Working Directory ",
                      control=control.CHOICE_BOX,
                      name='removeWorkingDir',
                      type=dtype.BOOL,
                      values=["On", "Off"],
                      comment="Deletes the contents of the Working "
                              "Directory after running.\n\nThis saves disk "
                              "space, but any additional preprocessing or "
                              "analysis will have to be completely re-run.")

        self.page.add(label="Run Logging ",
                      control=control.CHOICE_BOX,
                      name="run_logging",
                      type=dtype.BOOL,
                      values=["On", "Off"],
                      comment="Whether to write log details of the pipeline. "
                              "run to the logging files.")

        self.page.add(label="Regenerate Outputs ",
                      control=control.CHOICE_BOX,
                      name='reGenerateOutputs',
                      type=dtype.BOOL,
                      values=["On", "Off"],
                      comment="Uses the contents of the Working Directory "
                              "to regenerate all outputs and their "
                              "symbolic links.\n\nRequires an intact "
                              "Working Directory from a previous CPAC run.")

        self.page.add(label="Create Symbolic Links ",
                      control=control.CHOICE_BOX,
                      name='runSymbolicLinks',
                      type=dtype.LSTR,
                      comment="Create a user-friendly, well organized "
                              "version of the output directory.",
                      values=["On", "Off"])

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter


class WorkflowConfig(wx.ScrolledWindow):

    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)
        self.counter = counter

        self.page = GenericClass(self, "Preprocessing Workflow Selection")

        self.page.add(label="Gather Anatomical Data ",
                      control=control.CHOICE_BOX,
                      name='runAnatomicalDataGathering',
                      type=dtype.LSTR,
                      comment="Loads anatomical images for processing by CPAC.\n\nMust be enabled to run preprocessing and analyses.",
                      values=["On", "Off"])

        self.page.add(label="Gather Functional Data ",
                      control=control.CHOICE_BOX,
                      name='runFunctionalDataGathering',
                      type=dtype.LSTR,
                      comment="Loads functional images for processing by CPAC.\n\nMust be enabled to run preprocessing and analyses.",
                      values=["On", "Off"])

        self.page.add(label="Run Anatomical Preprocessing ",
                      control=control.CHOICE_BOX,
                      name='runAnatomicalPreprocessing',
                      type=dtype.LSTR,
                      comment="Runs the anatomical preprocessing workflow.\n\nMust be enabled to run any subsequent processing or analysis workflows.",
                      values=["On", "Off"])

        self.page.add(label="Run Functional Preprocessing ",
                      control=control.CHOICE_BOX,
                      name='runFunctionalPreprocessing',
                      type=dtype.LSTR,
                      comment="Runs the functional preprocessing workflow.\n\nMust be enabled to run any subsequent processing or analysis workflows.",
                      values=["On", "Off"])

        self.page.add(label="Functional Masking ",
                      control=control.CHOICE_BOX,
                      name='functionalMasking',
                      type=dtype.LSTR,
                      comment="Choose which tool to be used in functional masking - AFNI 3dAutoMask or FSL BET.",
                      values=["3dAutoMask", "BET", "3dAutoMask & BET"])

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter

