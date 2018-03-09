#!/usr/bin/python

# TODO: use Python Click to make this nice, or if not, use argparse

# TODO: probably need Click to nest the analysis presets and then THEIR
# TODO: required inputs


def main():

    import os
    import argparse

    from CPAC.utils import create_group_analysis_files

    parser = argparse.ArgumentParser()

    parser.add_argument("cpac_outputs", type=str,
                        help="the path to the CPAC pipeline output directory "
                             "containing individual-level analysis output "
                             "files for each participant")
    parser.add_argument("analysis_preset", type=str,
                        help="the type of group-level analysis to run using "
                             "FSL FLAME\n\nOptions:\n"
                             "single_grp: Single Group Average\n"
                             "single_grp_cov: Single Group Average w/ "
                             "Covariate\n"
                             "unpaired_two: Unpaired Two-Group Difference "
                             "(two-sample unpaired T-test)")
    parser.add_argument("z_thresh", type=str,
                        help="the z-threshold")
    parser.add_argument("p_thresh", type=str,
                        help="the p-threshold (cluster significance)")
    parser.add_argument("model_name", type=str,
                        help="name for the model")
    parser.add_argument("derivatives", type=str,
                        help="list of derivative names separated by spaces")
    parser.add_argument("--include", type=str, default=None,
                        help="the path to the group-level analysis "
                             "participant list text file")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="the output directory")
    parser.add_argument("--pheno_file", type=str, default=None,
                        help="the additional covariate for single-group "
                             "average")
    parser.add_argument("--pheno_sub_label", type=str, default=None,
                        help="the additional covariate for single-group "
                             "average")
    parser.add_argument("--covariate", type=str, default=None,
                        help="the additional covariate for single-group "
                             "average")

    args = parser.parse_args()

    if not args.output_dir:
        output_dir = os.getcwd()
    else:
        output_dir = args.output_dir

    if not args.include:
        include = [x for x in os.listdir(args.cpac_outputs) if os.path.isdir(x)]
    else:
        include = args.include

    derivatives_list = args.derivatives.split(" ")

    create_group_analysis_files.run(include, derivatives_list, args.z_thresh,
                                    args.p_thresh, args.analysis_preset,
                                    args.pheno_file, args.pheno_sub_label,
                                    output_dir, args.model_name,
                                    args.covariate)


if __name__ == "__main__":
    main()
