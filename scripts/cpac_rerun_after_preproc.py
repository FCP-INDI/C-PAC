#!/usr/bin/python

working_dir_keep = [
    '_scan_',
    'anat_gather_',
    'anat_mni_ants_register_',
    'anat_symmetric_mni_ants_register_',
    'anat_mni_fnirt_register_',
    'anat_preproc_',
    'anat_symmetric_mni_fnirt_register_',
    'd3.js',
    'edit_func_',
    'func_gather_',
    'func_preproc_automask_',
    'func_to_anat_bbreg_',
    'func_to_anat_FLIRT_',
    'graph.dot',
    'graph.dot.png',
    'graph.json',
    'graph1.json',
    'graph_detailed.dot',
    'graph_detailed.dot.png',
    'index.html',
    'log_anat_mni_fnirt_register_',
    'log_anat_preproc_',
    'log_anat_symmetric_mni_fnirt_register_',
    'log_fristons_parameter_model_',
    'log_func_preproc_automask_',
    'log_gen_motion_stats_',
    'log_motion_correct_to_standard_smooth_',
    'log_seg_preproc_',
    'seg_preproc_'
]


def main():
    """Clear a CPAC working directory of all node folders related to pipeline
    steps which occur after pre-processing.

    This script allows easy clearing of the working directory in order to
    re-run a pipeline without re-running registration and segmentation, which
    are time-consuming.
    """

    import os
    import shutil
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("working_directory", type=str,
                        help="file path to the working directory of the CPAC "
                             "run")

    args = parser.parse_args()

    work_dir = os.path.abspath(args.working_directory)

    for sub_id in os.listdir(work_dir):
        try:
            folders = os.listdir(os.path.join(work_dir, sub_id))
        except OSError:
            continue
        for dirname in folders:
            keep = False
            for substring in working_dir_keep:
                if substring in dirname:
                    keep = True

            if not keep:
                try:
                    shutil.rmtree(os.path.join(work_dir, dirname))
                    print ('Clearing {0} from the working '
                           'directory..'.format(dirname))
                except OSError:
                    pass

    print('Finished resetting working directory')


if __name__ == "__main__":
    main()
