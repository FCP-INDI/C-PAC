#!/usr/bin/python


def main():

    import os
    import sys
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("pipeline_config", type=str,
                        help="the path to the pipeline configuration YAML "
                             "file")
    parser.add_argument("data_config", type=str,
                        help="the path to the data configuration YAML file "
                             "(the participant list)")

    parser.add_argument("--nipype_install", type=str, default=None,
                        help="the build directory of a custom Nipype "
                             "installation you wish to use for this run")
    parser.add_argument("--cpac_install", type=str, default=None,
                        help="the build directory of a custom CPAC "
                             "installation you wish to use for this run")

    args = parser.parse_args()

    # put custom install directories at the beginning of sys.path
    if args.nipype_install:
        sys.path.insert(0, args.nipype_install)
    if args.cpac_install:
        sys.path.insert(0, args.cpac_install)

    if not os.path.exists(args.pipeline_config):
        err = "\n[!] The pipeline configuration file you provided " \
              "does not exist:\n{0}\n".format(args.pipeline_config)
        raise Exception(err)

    if not os.path.exists(args.data_config):
        err = "\n[!] The data configuration file you provided " \
              "does not exist:\n{0}\n".format(args.data_config)
        raise Exception(err)

    import CPAC
    CPAC.pipeline.cpac_runner.run(args.pipeline_config, args.data_config)


if __name__ == "__main__":
    main()
