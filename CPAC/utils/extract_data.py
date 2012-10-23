import sys
import os
import glob
import string
import logging

logging.basicConfig(filename='extract_data_logs.log', filemode='w', level=logging.DEBUG,\
                    format="%(levelname)s %(asctime)s %(lineno)d %(message)s")


def extract_data(c, param_map):
    """
    Method to generate a CPAC input subject list
    python file. The method extracts anatomical
    and functional data for each site( if multiple site)
    and/or scan parameters for each site and put it into
    a data structure read by python

    Example:
    subjects_list =[
       {
        'subject_id' : '0050386',
        'unique_id' : 'session_1',
        'anat': '/Users/home/data/NYU/0050386/session_1/anat_1/anat.nii.gz',
        'rest':{
            'rest_1_rest' : '/Users/home/data/NYU/0050386/session_1/rest_1/rest.nii.gz',
            'rest_2_rest' : '/Users/home/data/NYU/0050386/session_1/rest_2/rest.nii.gz',
            }
        'scan_parameters':{
            'tr': '2',
            'acquisition': 'alt+z2',
            'reference': '17',
            'first_tr': '',
            'last_tr': '',
            }
        },
    ]

    or

    subjects_list =[
       {
        'subject_id' : '0050386',
        'unique_id' : 'session_1',
        'anat': '/Users/home/data/NYU/0050386/session_1/anat_1/anat.nii.gz',
        'rest':{
            'rest_1_rest' : '/Users/home/data/NYU/0050386/session_1/rest_1/rest.nii.gz',
            'rest_2_rest' : '/Users/home/data/NYU/0050386/session_1/rest_2/rest.nii.gz',
            }
          },
    ]

    """

    #method to read each line of the file into list
    #returns list
    def get_list(arg):
        if isinstance(arg, list):
            ret_list = arg
        else:
            ret_list = [fline.rstrip('\r\n') for fline in open(arg, 'r').readlines()]

        return ret_list

    exclusion_list = []
    if c.exclusionSubjectList is not None:
        exclusion_list = get_list(c.exclusionSubjectList)

    subject_list = []
    if c.subjectList is not None:
        subject_list = get_list(c.subjectList)

    #check if Template is correct
    def checkTemplate(template):

        if template.count('%s') != 2:
            logging.exception("Please provide '%s' in the template" \
                            "where your site and subjects are present"\
                            "Please see examples")

        filename, ext = os.path.splitext(os.path.basename(template))
        ext = os.path.splitext(filename)[1] + ext

        if ext not in [".nii", ".nii.gz"]:
            logging.exception("Invalid file name", os.path.basename(template))

    def get_site_list(path):
        base, relative = path.split('%s')
        sites = os.listdir(base)
        return sites

    def create_site_subject_mapping(base, relative):

        #mapping between site and subject
        site_subject_map = {}
        base_path_list = []

        if c.siteList is not None:
            site_list = get_list(c.siteList)
        else:
            site_list = get_site_list(base)

        for site in site_list:
            paths = glob.glob(string.replace(base, '%s', site))
            base_path_list.extend(paths)
            for path in paths:
                for sub in os.listdir(path):
                    #check if subject is present in subject_list
                    if subject_list:
                        if sub in subject_list:
                            site_subject_map[sub] = site
                    elif sub not in exclusion_list:
                        if sub not in '.DS_Store':
                            site_subject_map[sub] = site

        return base_path_list, site_subject_map

    #method to split the input template path
    #into base, path before subject directory
    #and relative, path after subject directory
    def getPath(template):

        checkTemplate(template)
        base, relative = template.rsplit("%s", 1)
        base, subject_map = create_site_subject_mapping(base, relative)
        base.sort()
        relative = relative.lstrip("/")
        return base, relative, subject_map

    #get anatomical base path and anatomical relative path
    anat_base, anat_relative = getPath(c.anatomicalTemplate)[:2]

    #get functional base path, functional relative path and site-subject map
    func_base, func_relative, subject_map = getPath(c.functionalTemplate)

    if not anat_base:
        logging.exception("No such file or directory %s", anat_base)
        logging.exception("Anatomical Data template incorrect")

    if not func_base:
        logging.exception("No such file or directory %s", func_base)
        logging.exception("Functional Data template incorrect")

    if len(anat_base) != len(func_base):
        logging.exception("Some sites are missing, Please check your template"\
              , anat_base, "!=", func_base)
        logging.exception(" Base length Unequal. Some sites are missing."\
                           "extract_data doesn't script support this.Please" \
                           "Provide your own subjects_list file")

    #calculate the length of relative paths(path after subject directory)
    func_relative_len = len(func_relative.split('/'))
    anat_relative_len = len(anat_relative.split('/'))

    def check_for_sessions(relative_path, path_length):
        """
        Method to check if there are sessions present

        """
        #default
        session_present = False
        session_path = 'session_1'

        #session present if path_length is equal to 3
        if path_length == 3:
            relative_path_list = relative_path.split('/')
            session_path = relative_path_list[0]
            relative_path = string.join(relative_path_list[1:], "/")
            session_present = True
        elif path_length > 3:
            logging.exception("extract_data script currently doesn't support this directory structure."\
                             "Please provide the subjects_list file to run CPAC."\
                             "For more information refer to manual")
        return session_present, session_path, relative_path

    func_session_present, func_session_path, func_relative = \
        check_for_sessions(func_relative, func_relative_len)

    anat_session_present, anat_session_path, anat_relative = \
        check_for_sessions(anat_relative, anat_relative_len)

    f = open("CPAC_subject_list.py", 'wb')

    def fetch_path(i, anat_sub, func_sub, session_id):
        """
        Method to extract anatomical and functional
        path for a session and print to file

        Parameters
        ----------
        i : int
            index of site
        anat_sub : string
            string containing subject/ concatenated
            subject-session path for anatomical file
        func_sub: string
            string containing subject/ concatenated
            subject-session path for functional file
        session_id: string
            session

        Raises
        ------
        Exception
        """

        try:

            def print_begin_of_file(sub, session_id):
                print >> f, "{"
                print >> f, "    'subject_id': '" + sub + "',"
                print >> f, "    'unique_id': '" + session_id + "',"

            def print_end_of_file(sub):
                if param_map is not None:
                    try:
                        logging.debug("site for sub %s -> %s", sub, subject_map.get(sub))
                        logging.debug("scan parameters for the above site %s", param_map.get(subject_map.get(sub)))
                        print >> f, "    'scan_parameters':{"
                        print >> f, "        'tr': '" + param_map.get(subject_map.get(sub))[4] + "',"
                        print >> f, "        'acquisition': '" + param_map.get(subject_map.get(sub))[0] + "',"
                        print >> f, "        'reference': '" + param_map.get(subject_map.get(sub))[3] + "'" + ","
                        print >> f, "        'first_tr': '" + param_map.get(subject_map.get(sub))[1] +  "'" + ","
                        print >> f, "        'last_tr': '" + param_map.get(subject_map.get(sub))[2] + "'" + ","
                        print >> f, "        }"
                    except:
                        logging.exception(" No Parameter values for the %s site is defined in the scan"\
                                        " parameters csv file", subject_map.get(sub))

                print >> f, "},"

            #get anatomical file
            anat_base_path = os.path.join(anat_base[i], anat_sub)
            func_base_path = os.path.join(func_base[i], func_sub)

            anat = None
            func = None

            anat = glob.glob(os.path.join(anat_base_path, anat_relative))
            func = glob.glob(os.path.join(func_base_path, func_relative))

            if anat and func:
                print_begin_of_file(anat_sub.split("/")[0], session_id)
                print >> f, "    'anat': '" + anat[0] + "',"
                print >>f, "    'rest':{"

                #iterate for each rest session
                for iter in func:
                    #get scan_id
                    iterable = os.path.splitext(os.path.splitext(iter.replace(func_base_path,'').lstrip("/"))[0])[0]
                    iterable = iterable.replace("/", "_")
                    print>>f,  "      '" + iterable + "': '" + iter + "',"
                print >> f, "      },"
                print_end_of_file(anat_sub.split("/")[0])
            else:
                logging.debug("skipping subject %s", anat_sub.split("/")[0])

        except Exception:
            logging.exception("Exception while fetching anatomical and functional paths")

    def walk(index, sub):
        """
        Method which walks across each subject
        path in the data site path

        Parameters
        ----------
        index : int
            index of site
        sub : string
            subject_id

        Raises
        ------
        Exception
        """
        try:

            if func_session_present:
                #if there are sessions
                if "*" in func_session_path:
                    session_list = glob.glob(os.path.join(func_base[index],os.path.join(sub, func_session_path)))
                    if session_list:
                        for session in session_list:
                            session_id = os.path.basename(session)
                            if anat_session_present:
                                if func_session_path == anat_session_path:
                                    fetch_path(index, os.path.join(sub,session_id), os.path.join(sub,session_id), session_id)
                                else:
                                    fetch_path(index, os.path.join(sub, anat_session_path), os.path.join(sub, session_id), session_id)
                            else:
                                fetch_path(index, sub, os.path.join(sub, session_id), session_id)
                    else:
                        logging.debug("Skipping subject %s", sub)
                else:
                    session_id = func_session_path
                    fetch_path(index, os.path.join(sub, anat_session_path), os.path.join(sub, func_session_path), session_id) 
            else:
                logging.debug("No sessions")
                session_id = ''
                fetch_path(index, sub, sub, session_id)

        except Exception:
            logging.exception(Exception.message)
        except:
            logging.exception("Please make sessions are consistent across all subjects")

    try:
        print >>f, "subjects_list = ["
        for i in range(len(anat_base)):
            for sub in os.listdir(anat_base[i]):
                #check if subject is present in subject_list
                if subject_list:
                    if sub in subject_list:
                        logging.debug("extracting data for subject: %s", sub)
                        walk(i, sub)
                #check that subject is not in exclusion list
                elif sub not in exclusion_list and sub not in '.DS_Store':
                    logging.debug("extracting data for subject: %s", sub)
                    walk(i, sub)

        print >> f, "]"

        name = os.path.join(os.getcwd(), 'CPAC_subject_list.py')
        print "Extraction Successfully Completed...Input Subjects_list for CPAC - %s" %name
    except Exception:
        logging.exception(Exception.message)
    finally:
        f.close()


def generate_suplimentary_files():
    """
    Method to generate phenotypic template file
    and subject list for group analysis
    """
    from sets import Set
    import csv

    c = __import__('CPAC_subject_list')

    subject_scan_set = Set()
    subject_set = Set()
    scan_set = Set()
    data_list = []

    for sub in c.subjects_list:
        subject_id = sub['subject_id'] + "_" + sub['unique_id']
        for scan in sub['rest'].keys():
            subject_scan_set.add((subject_id, scan))
            subject_set.add(subject_id)
            scan_set.add(scan)

    for item in subject_scan_set:
        list1 = []
        list1.append(item[0] + "/" + item[1])
        for val in subject_set:
            if val in item:
                list1.append(1)
            else:
                list1.append(0)

        for val in scan_set:
            if val in item:
                list1.append(1)
            else:
                list1.append(0)

        data_list.append(list1)

    #prepare data for phenotypic file
    if len(scan_set) > 1:
        list1 = ['subject_id/scan']
        list1.extend(list(subject_set))
        list1.extend(list(scan_set))

    file_name = os.path.join(os.getcwd(),'phenotypic_template.csv')
    f = open(file_name, 'wb')
    writer = csv.writer(f)

    if len(scan_set) > 1:
        writer.writerow(list1)
        writer.writerows(data_list)
    else:
        writer.writerow(['subject_id'])
        for sub in subject_set:
            writer.writerow([sub])

    f.close()

    print "Template Phenotypic file for group analysis - %s" %file_name

    file_name = os.path.join(os.getcwd(), "subject_list_group_analysis.txt")
    f = open(file_name, 'w')

    for sub in subject_set:
        print >> f, sub

    print "Subject list required later for group analysis - %s" %file_name
    f.close()


def read_csv(csv_input):
    """
    Method to read csv file
     'Acquisition'
     'Reference'
     'Site'
     'TR (seconds)'

    """
    import csv
    from collections import defaultdict
    try:
        reader = csv.DictReader(open(csv_input, "U"))

        dict_labels = defaultdict(list)
        for line in reader:
            csv_dict = dict((k.lower(), v) for k, v in line.iteritems())
            dict_labels[csv_dict.get('site')] = [csv_dict[key] for key in sorted(csv_dict.keys()) \
                                                 if key != 'site' and key != 'scan']

        if len(dict_labels.keys()) < 1:
            logging.exception("Scan Parameters File is either empty"\
                              "or missing header")

        return dict_labels

    except IOError:
        logging.exception("Error reading the csv file %s", csv_input)
    except:
        logging.exception("Error reading scan parameters csv. Make sure you are using the correct template")


def run(data_config):
    """
    Run method takes data_config
    file as the input argument
    """

    path, fname = os.path.split(os.path.realpath(data_config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])
    print "For any errors or messages check the log file - %s"\
           % os.path.join(os.getcwd(), 'extract_data_logs.log')

    if c.scanParametersCSV is not None:
        s_param_map = read_csv(c.scanParametersCSV)
    else:
        logging.debug("no scan parameters csv included"\
              "make sure you turn off slice timing correction option"\
              "in CPAC configuration")
        s_param_map = None

    extract_data(c, s_param_map)
    generate_suplimentary_files()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: python extract_data.py data_config.py"
        sys.exit()
    else:
        run(sys.argv[1])