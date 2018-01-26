#!/usr/bin/env python


def print_data_config_info(data_config_yml):
    # TODO: implement this!
    pass


def main():

    import os
    import CPAC
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--data_settings_file", type=str, default=None,
                        help="the path to the data settings YML file, which "
                             "contains information like the data format, "
                             "BIDS base directory, or anatomical and "
                             "functional file templates")
    parser.add_argument("--generate_template",
                        action='store_true', default=False,
                        help="create a blank template of the data settings "
                             "file - configure your settings and provide "
                             "this file to this script to generate your "
                             "data configuration file")
    '''
    parser.add_argument("--data_config_info", type=str, default=None,
                        help="the path to a data configuration YML file "
                             "you want information about")
    '''

    args = parser.parse_args()

    if not args.data_settings_file and not args.generate_template: # and \
            #not args.data_config_info:
        print "No inputs provided. Use the -h flag for instructions.\n"

    if args.data_settings_file and not args.generate_template: # and \
            #not args.data_config_info:

        input_file = os.path.abspath(args.data_settings_file)
        if os.path.exists(input_file):
            if not input_file.endswith((".yml", ".yaml", ".YML", ".YAML")):
                err = "\n[!] Data settings file must be a YAML " \
                      "(.yml/.yaml) file.\n"
                print err
            else:
                # create the data config!
                CPAC.utils.build_data_config.run(input_file)
        else:
            err = "\n[!] Data settings file path cannot be found:\n" \
                  "{0}\n".format(input_file)
            print err

    elif args.generate_template and not args.data_settings_file: #and \
            #not args.data_config_info:

        import shutil
        import pkg_resources as p
        settings_template = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                             "configs",
                                             "data_settings_template.yml"))

        settings_file = os.path.join(os.getcwd(), "data_settings.yml")

        try:
            if os.path.exists(settings_file):
                settings_file = os.path.join(os.getcwd(),
                                             "data_settings_1.yml")
                while os.path.exists(settings_file):
                    idx = int(os.path.basename(settings_file).split("_")[2].replace(".yml", ""))
                    settings_file = os.path.join(os.getcwd(),
                                                 "data_settings_{0}.yml".format(idx+1))
            shutil.copy(settings_template, settings_file)
        except:
            err = "\n[!] Could not write the data settings file template " \
                  "to the current directory.\n"
            raise Exception(err)

        print "\nGenerated a default data settings YML file for editing:\n" \
              "{0}\n\nThis file can be completed and entered into the " \
              "cpac_data_config_setup.py script with the " \
              "--data_settings_file flag.\nAdditionally, it can also be " \
              "loaded into the CPAC data configuration file builder UI " \
              "using the 'Load Preset' button.\n".format(settings_file)

    else:
        print "Too many arguments. Only one option is accepted at a time.\n"

    '''
    elif args.data_config_info and not args.data_settings_file and \
            not args.generate_template:

        input_file = os.path.abspath(args.data_config_info)
        if os.path.exists(input_file):
            if not input_file.endswith((".yml", ".yaml", ".YML", ".YAML")):
                err = "\n[!] Data configuration file must be a YAML " \
                      "(.yml/.yaml) file.\n"
                print err
            else:
                # get data config info!!
                print_data_config_info(input_file)

        else:
            err = "\n[!] Data configuration file path cannot be found:\n" \
                  "{0}\n".format(input_file)
            print err

    else:
        print "Too many arguments. Only one option is accepted at a time.\n"
    '''


if __name__ == "__main__":
    main()
