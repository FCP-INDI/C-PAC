#!/usr/bin/env python


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
                        help="Write a summary report in PDF "
                        "format.")

    args = parser.parse_args()

    #TODO: more modular/non-redundant way to do this?
    settings_text = \
        ""

    if not args.data_settings_file and not args.generate_template:
        print "No inputs provided. Use the -h flag for instructions.\n"

    if args.data_settings_file:
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

    if args.generate_template:
        settings_file = os.path.join(os.getcwd(), "cpac_data_settings.yml")

        try:
            with open(settings_file, "wt") as f:
                f.write(settings_text)
        except:
            err = "\n[!] Could not write the data settings file template " \
                  "to the current directory.\n"
            raise Exception(err)

        print "\nGenerated a default data settings YML file for editing:\n" \
              "{0}\n\nThis file can be completed and entered into the " \
              "cpac_data_config_setup.py script with the " \
              "--data_settings_file flag.\n".format(settings_file)


if __name__ == "__main__":
    main()
