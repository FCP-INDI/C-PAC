"""
Tests for surface configuration

FROM: abcd-options
pipeline_setup:
  pipeline_name: ABCD-blip
  output_directory:
    quality_control:
      generate_xcpqc_files: On
functional_preproc:
  distortion_correction:
    using:
      - PhaseDiff
      - Blip
"""