import sys
import os
import glob
import string
import logging
import yaml

#logging.basicConfig(filename=os.path.join(os.getcwd(), 'extract_data_logs.log'), filemode='w', level=logging.DEBUG,\
#                    format="%(levelname)s %(asctime)s %(lineno)d %(message)s")


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
            msg = "Please provide '%s' in the template" \
                  "where your site and subjects are present"\
                  "Please see examples"
            logging.exception(msg)
            raise Exception(msg)

        filename, ext = os.path.splitext(os.path.basename(template))
        ext = os.path.splitext(filename)[1] + ext

        if ext not in [".nii", ".nii.gz"]:
            msg = "Invalid file name", os.path.basename(template)
            logging.exception(msg)
            raise Exception(msg)

    def get_site_list(path):
        base, relative = path.split('%s')
        sites = os.listdir(base)
        return sites
    
    def check_length(scan_name, file_name):
        
        if len(file_name) > 30:
            msg = "filename- %s is too long."\
                   "It should not be more than 30 characters."%(file_name)
            logging.exception(msg)
            raise Exception(msg)
        
        if len(scan_name) - len(os.path.splitext(os.path.splitext(file_name)[0])[0])>= 40:
            msg = "scan name %s is too long."\
                  "It should not be more than 20 characters"\
                  %(scan_name.replace("_"+os.path.splitext(os.path.splitext(file_name)[0])[0], ''))
            logging.exception(msg)
            raise Exception(msg)
        

        

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
                        if sub in subject_list and sub not in exclusion_list:
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
        msg = "Anatomical Data template incorrect. No such file or directory %s", anat_base
        logging.exception(msg)
        raise Exception(msg)

    if not func_base:
        msg = "Functional Data template incorrect. No such file or directory %s, func_base"
        logging.exception(msg)
        raise Exception(msg)
        
    if len(anat_base) != len(func_base):
        msg1 = "Some sites are missing, Please check your template"\
              , anat_base, "!=", func_base
        logging.exception(msg1)
        
        msg2 = " Base length Unequal. Some sites are missing."\
               "extract_data doesn't script support this.Please" \
               "Provide your own subjects_list file"
        logging.exception(msg2)
        raise Exception(msg2)

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
            msg = "extract_data script currently doesn't support this directory structure."\
                  "Please provide the subjects_list file to run CPAC."\
                  "For more information refer to manual"
            logging.exception(msg)
            raise Exception(msg)
        return session_present, session_path, relative_path

    func_session_present, func_session_path, func_relative = \
        check_for_sessions(func_relative, func_relative_len)

    anat_session_present, anat_session_path, anat_relative = \
        check_for_sessions(anat_relative, anat_relative_len)

    f = open(os.path.join(c.outputSubjectListLocation, "CPAC_subject_list_%s.yml" % c.subjectListName[0]), 'wb')



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
                print >> f, "-"
                print >> f, "    subject_id: '" + sub + "'"
                print >> f, "    unique_id: '" + session_id + "'" 

            def print_end_of_file(sub):
                if param_map is not None:
                    try:
                        logging.debug("site for sub %s -> %s" %(sub, subject_map.get(sub)))
                        logging.debug("scan parameters for the above site %s"%param_map.get(subject_map.get(sub)))
                        print >> f, "    scan_parameters:"
                        print >> f, "        tr: '" + param_map.get(subject_map.get(sub))[4] + "'" 
                        print >> f, "        acquisition: '" + param_map.get(subject_map.get(sub))[0] + "'" 
                        print >> f, "        reference: '" + param_map.get(subject_map.get(sub))[3] + "'"
                        print >> f, "        first_tr: '" + param_map.get(subject_map.get(sub))[1] + "'"
                        print >> f, "        last_tr: '" + param_map.get(subject_map.get(sub))[2] + "'"
                    except:
                        msg = " No Parameter values for the %s site is defined in the scan"\
                              " parameters csv file" %subject_map.get(sub)
                        raise ValueError(msg)

            #get anatomical file
            anat_base_path = os.path.join(anat_base[i], anat_sub)
            func_base_path = os.path.join(func_base[i], func_sub)

            anat = None
            func = None

            anat = glob.glob(os.path.join(anat_base_path, anat_relative))
            func = glob.glob(os.path.join(func_base_path, func_relative))

            if anat and func:
                print_begin_of_file(anat_sub.split("/")[0], session_id)
                print >> f, "    anat: '" + os.path.realpath(anat[0]) + "'"
                print >> f, "    rest: "

                #iterate for each rest session
                for iter in func:
                    #get scan_id
                    iterable = os.path.splitext(os.path.splitext(iter.replace(func_base_path, '').lstrip("/"))[0])[0]
                    iterable = iterable.replace("/", "_")
                    check_length(iterable, os.path.basename(os.path.realpath(iter)))
                    print>>f, "      " + iterable + ": '" + os.path.realpath(iter) + "'"
                
                print_end_of_file(anat_sub.split("/")[0])
                
            else:
                logging.debug("skipping subject %s"%anat_sub.split("/")[0])
        
        except ValueError:

            logging.exception(ValueError.message)
            raise

        except Exception, e:

            err_msg = 'Exception while felching anatomical and functional ' \
                      'paths: \n' + str(e)

            logging.exception(err_msg)
            raise Exception(err_msg)



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
                    session_list = glob.glob(os.path.join(func_base[index], os.path.join(sub, func_session_path)))
                else:
                    session_list = [func_session_path]

                if session_list:
                    for session in session_list:
                        session_id = os.path.basename(session)
                        if anat_session_present:
                            if func_session_path == anat_session_path:
                                fetch_path(index, os.path.join(sub, session_id), os.path.join(sub, session_id), session_id)
                            else:
                                fetch_path(index, os.path.join(sub, anat_session_path), os.path.join(sub, session_id), session_id)
                        else:
                            fetch_path(index, sub, os.path.join(sub, session_id), session_id)
                else:
                    logging.debug("Skipping subject %s", sub)

            else:
                logging.debug("No sessions")
                session_id = ''
                fetch_path(index, sub, sub, session_id)

        except Exception:

            logging.exception(Exception.message)
            raise

        except:

            err_msg = 'Please make sessions are consistent across all ' \
                      'subjects.\n\n'

            logging.exception(err_msg)
            raise Exception(err_msg)


    try:
        for i in range(len(anat_base)):
            for sub in os.listdir(anat_base[i]):
                #check if subject is present in subject_list
                if subject_list:
                    if sub in subject_list and sub not in exclusion_list:
                        logging.debug("extracting data for subject: %s", sub)
                        walk(i, sub)
                #check that subject is not in exclusion list
                elif sub not in exclusion_list and sub not in '.DS_Store':
                    logging.debug("extracting data for subject: %s", sub)
                    walk(i, sub)

        
        name = os.path.join(c.outputSubjectListLocation, 'CPAC_subject_list.yml')
        print "Extraction Successfully Completed...Input Subjects_list for CPAC - %s" % name

    except Exception:

        logging.exception(Exception.message)
        raise

    finally:

        f.close()




def generate_supplementary_files(output_path, subject_list_name):

    """
    Method to generate phenotypic template file
    and subject list for group analysis
    """

    from sets import Set
    import csv

    subject_list_name = subject_list_name[0]

    try:
        subjects_list = yaml.load(open(os.path.join(output_path, 'CPAC_' \
                'subject_list_%s.yml' % subject_list_name), 'r'))
    except:
        print 'Subject list couldn\'t be read!'
        print 'path: ', os.path.join(output_path, 'CPAC_subject_list_%s.yml' \
                % subject_list_name)
        raise Exception

    subject_scan_set = Set()
    subID_set = Set()
    session_set = Set()
    subject_set = Set()
    scan_set = Set()
    data_list = []

    try:
        for sub in subjects_list:
            
            if sub['unique_id']:
                subject_id = sub['subject_id'] + "_" + sub['unique_id']
            else:
                subject_id = sub['subject_id']
                
            for scan in sub['rest'].keys():
                subject_scan_set.add((subject_id, scan))
                subID_set.add(sub['subject_id'])
                session_set.add(sub['unique_id'])
                subject_set.add(subject_id)
                scan_set.add(scan)
    except TypeError as e:
        print 'Subject list could not be populated!'
        print 'This is most likely due to a mis-formatting in your '\
              'inclusion and/or exclusion subjects txt file or your '\
              'anatomical and/or functional path templates.'
        print 'Error: %s' % e
        err_str = 'Check formatting of your anatomical/functional path '\
                  'templates and inclusion/exclusion subjects text files'
        raise TypeError(err_str)

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



    # generate the phenotypic file templates for group analysis

    file_name = os.path.join(output_path, 'phenotypic_template_%s.csv' \
            % subject_list_name)

    try:
        f = open(file_name, 'wb')
    except:
        print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
        print file_name, '\n\n'
        print 'Make sure you have write access? Then come back. Don\'t ' \
                'worry.. I\'ll wait.\n\n'
        raise IOError

    writer = csv.writer(f)

    writer.writerow(['subject_id', 'EV1', '..'])
    for sub in sorted(subID_set):
        writer.writerow([sub, ''])

    f.close()

    print "Template Phenotypic file for group analysis - %s" % file_name


    # generate the phenotypic file templates for repeated measures
    if (len(session_set) > 1) and (len(scan_set) > 1):

        file_name = os.path.join(output_path, 'phenotypic_template_repeated' \
                '_measures_mult_sessions_and_scans_%s.csv' \
                % subject_list_name)

        try:
            f = open(file_name, 'wb')
        except:
            print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
            print file_name, '\n\n'
            print 'Make sure you have write access? Then come back. Don\'t ' \
                    'worry.. I\'ll wait.\n\n'
            raise IOError

        writer = csv.writer(f)

        writer.writerow(['subject_id', 'EV1', '..'])

        for session in sorted(session_set):
            for scan in sorted(scan_set):
                for sub in sorted(subID_set):
                    writer.writerow([sub + '_' + scan + '_' + session, ''])

        f.close()



    if (len(session_set) > 1):

        file_name = os.path.join(output_path, 'phenotypic_template_repeated' \
                '_measures_multiple_sessions_%s.csv' % subject_list_name)

        try:
            f = open(file_name, 'wb')
        except:
            print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
            print file_name, '\n\n'
            print 'Make sure you have write access? Then come back. Don\'t ' \
                    'worry.. I\'ll wait.\n\n'
            raise IOError

        writer = csv.writer(f)

        writer.writerow(['subject_id', 'EV1', '..'])

        for session in sorted(session_set):
            for sub in sorted(subID_set):
                writer.writerow([sub + '_' + session, ''])

        f.close()



    if (len(scan_set) > 1):

        file_name = os.path.join(output_path, 'phenotypic_template_repeated' \
                '_measures_multiple_scans_%s.csv' % subject_list_name)

        try:
            f = open(file_name, 'wb')
        except:
            print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
            print file_name, '\n\n'
            print 'Make sure you have write access? Then come back. Don\'t ' \
                    'worry.. I\'ll wait.\n\n'
            raise IOError

        writer = csv.writer(f)

        writer.writerow(['subject_id', 'EV1', '..'])

        for scan in sorted(scan_set):
            for sub in sorted(subID_set):
                writer.writerow([sub + '_' + scan, ''])

        f.close()




    # generate the group analysis subject lists

    file_name = os.path.join(output_path, 'subject_list_group_analysis' \
            '_%s.txt' % subject_list_name)

    try:
        f = open(file_name, 'w')
    except:
        print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
        print file_name, '\n\n'
        print 'Make sure you have write access? Then come back. Don\'t ' \
                'worry.. I\'ll wait.\n\n'
        raise IOError

    for sub in sorted(subID_set):
        print >> f, sub


    print "Subject list required later for group analysis - %s" % file_name
    f.close()


    # generate the group analysis subject lists for repeated measures
    if (len(session_set) > 1) and (len(scan_set) > 1):

        file_name = os.path.join(output_path, 'subject_list_group_analysis_' \
                'repeated_measures_mult_sessions_and_scans_%s.txt' \
                % subject_list_name)

        try:
            f = open(file_name, 'w')
        except:
            print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
            print file_name, '\n\n'
            print 'Make sure you have write access? Then come back. Don\'t ' \
                    'worry.. I\'ll wait.\n\n'
            raise IOError

        for session in sorted(session_set):
            for scan in sorted(scan_set):
                for sub in sorted(subID_set):
                    print >> f, sub + ',' + scan + ',' + session

        f.close()



    if (len(session_set) > 1):

        file_name = os.path.join(output_path, 'subject_list_group_analysis_' \
                'repeated_measures_multiple_sessions_%s.txt' \
                % subject_list_name)

        try:
            f = open(file_name, 'w')
        except:
            print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
            print file_name, '\n\n'
            print 'Make sure you have write access? Then come back. Don\'t ' \
                    'worry.. I\'ll wait.\n\n'
            raise IOError

        for session in sorted(session_set):
            for sub in sorted(subID_set):
                print >> f, sub + ',' + session

        f.close()



    if (len(scan_set) > 1):

        file_name = os.path.join(output_path, 'subject_list_group_analysis_' \
                'repeated_measures_multiple_scans_%s.txt' \
                % subject_list_name)

        try:
            f = open(file_name, 'w')
        except:
            print '\n\nCPAC says: I couldn\'t save this file to your drive:\n'
            print file_name, '\n\n'
            print 'Make sure you have write access? Then come back. Don\'t ' \
                    'worry.. I\'ll wait.\n\n'
            raise IOError

        for scan in sorted(scan_set):
            for sub in sorted(subID_set):
                print >> f, sub + ',' + scan

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
            msg ="Scan Parameters File is either empty"\
                 "or missing header"
            logging.exception(msg)
            raise Exception(msg)

        return dict_labels

    except IOError:
        msg = "Error reading the csv file %s", csv_input
        logging.exception(msg)
        raise Exception(msg)
    except:
        msg = "Error reading scan parameters csv. Make sure you are using the correct template"
        logging.exception(msg)
        raise Exception(msg)


"""
Class to set dictionary keys as map attributes
"""
class Configuration(object):
    def __init__(self, config_map):
        for key in config_map:
            if config_map[key] == 'None':
                config_map[key] = None
            setattr(self, key, config_map[key])
        

def run(data_config):
    """
    Run method takes data_config
    file as the input argument
    """
    root = logging.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    logging.basicConfig(filename=os.path.join(os.getcwd(), 'extract_data_logs.log'), filemode='w', level=logging.DEBUG,\
                    format="%(levelname)s %(asctime)s %(lineno)d %(message)s")

    print "For any errors or messages check the log file - %s"\
           % os.path.join(os.getcwd(), 'extract_data_logs.log')
    
    c = Configuration(yaml.load(open(os.path.realpath(data_config), 'r')))

    if c.scanParametersCSV is not None:
        s_param_map = read_csv(c.scanParametersCSV)
    else:
        logging.debug("no scan parameters csv included"\
              "make sure you turn off slice timing correction option"\
              "in CPAC configuration")
        s_param_map = None

    extract_data(c, s_param_map)
    generate_supplementary_files(c.outputSubjectListLocation, c.subjectListName)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: python extract_data.py data_config.yml"
        sys.exit()
    else:
        run(sys.argv[1])
