import argparse
import sys
import os
import glob
import string
           
def extract_data(c):
    """
    Method to create a python file
    containing subject list.
    The method extracts anatomical 
    and functional data for each site( if nultiple site)
    and put it into a data structure read by python
    
    Example:
    subjects_list =[
       {    
        'Subject_id' : '0050386',
        'Unique_id' : 'session_1',
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
    def get_list(fname):
        flines = open(fname, 'r').readlines()
        return [fline.rstrip('\r\n') for fline in flines]
 
    exclusion_list=[]
    
     
    if c.exclusionSubjectList is not None:
        exclusion_list= get_list(c.exclusionSubjectList)
    else:
        exclusion_list=[]


    if c.subjectList is not None:
        subject_list = get_list(c.subjectList)
    else:
        subject_list=[]

    #check if Template is correct
    def checkTemplate(template):
    
        if '%s' not in template:
            raise Exception("Please provide '%s' in the template" \
                            "where subjects are present")
    
        filename, ext = os.path.splitext(os.path.basename(template))
        ext = os.path.splitext(filename)[1] + ext
        
        if ext not in [".nii", ".nii.gz"]:
            raise Exception("Invalid file name",os.path.basename(template) )
        
        
    #method to split the input template path
    #into base, path before subject directory
    #and relative, path after subject directory
    def getPath(template):
        
        checkTemplate(template)
        base, relative = template.split("%s")
        base = glob.glob(base)
        base.sort()
        relative = relative.lstrip("/")
        return base, relative
    
    #get anatomical base path and anatomical relative path
    anat_base, anat_relative = getPath(c.anatomicalTemplate)  
    
    #get functional base path and functional relative path
    func_base, func_relative = getPath(c.functionalTemplate)
    
    if not anat_base:
        print "No such file or directory ", anat_base 
        raise Exception("Anotomical Data template incorrect")
    
    if not func_base:
        print "No such file or directory", func_base
        raise Exception ("Functional Data template incorrect")
    
    if len(anat_base) != len(func_base) :
        print "Some sites are missing, Please check your template" \
              ,anat_base,"!=", func_base 
        raise Exception (" Base length Unequal. Some sites are missing."\
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
            raise Exception("extract_data script currently doesn't support this directory structure."\
                             "Please provide the subjects_list file to run CPAC."\
                             "For more information refer to manual")    
            
        return session_present, session_path, relative_path
    
    if func_relative_len!= anat_relative_len:
        raise Exception(" extract_data script currently doesn't"\
                          "support different relative paths for"\
                          "Anatomical and functional files")
    
    func_session_present, func_session_path, func_relative = \
        check_for_sessions(func_relative, func_relative_len)
        
    anat_session_present, anat_session_path, anat_relative = \
        check_for_sessions(anat_relative, anat_relative_len)
    
    f = open("CPAC_subject_list.py", 'wb')
   
    
    def fetch_path(i, anat_sub, func_sub):
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
            
        Raises
        ------
        Exception
        """
        
        try:
            
            #get anatomical file
            anat_base_path = os.path.join(anat_base[i],anat_sub )

            if not os.path.exists(anat_base_path):
                print "path doesn't exist", anat_base_path
                raise Exception ("invalid Path. Please check anatomicalTemplate in the config file")
            
            anat = glob.glob(os.path.join(anat_base_path, anat_relative))[0]
            print >> f, "    'anat': '" + anat + "',"
             
            if not anat:
                print "Unable to find anatomical image at ",os.path.join(anat_base_path, anat_relative) 
                raise Exception("Anatomical Data Missing")
           
            #get functional file
            
            func_base_path = os.path.join(func_base[i], func_sub)
            
            if not os.path.exists(func_base_path):
                print "path doesn't exist", func_base_path
                raise Exception ("invalid Path. Please check functionalTemplate in the config file")
            
            func = glob.glob(os.path.join(func_base_path, func_relative))
           
            if not func:
                print "Unable to find teh functional image at", os.path.join(func_base_path, func_relative)
                raise Exception("Functional Data Missing")
           
            print >>f, "    'rest':{" 
            
            #iterate for each rest session 
            for iter in func :
                #get scan_id
                iterable = os.path.splitext(os.path.splitext(iter.replace(func_base_path,'').lstrip("/"))[0])[0]
                iterable = iterable.replace("/", "_")
                print>>f,  "      '"+iterable+"': '"+iter+"'," 
            print >> f, "      }"
            
    
        except Exception:
            raise
        
    
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
        
        def print_to_file(sub, session_id):
            print >> f, "{"
            print >> f, "    'Subject_id': '" + sub + "',"
            print >> f, "    'Unique_id': '" + session_id + "',"
        try:
            if func_session_present and anat_session_present:
                #if there are sessions
                if "*" in func_session_path:
                    session_list = glob.glob(os.path.join(func_base[index],os.path.join(sub, func_session_path)))
                    for session in session_list:
                        session_id= os.path.basename(session)
                        print_to_file(sub, session_id)
                        if func_session_path == anat_session_path:  
                            fetch_path(index, os.path.join(sub,session_id), os.path.join(sub,session_id))
                        else:
                            fetch_path(index, os.path.join(sub, anat_session_path), os.path.join(sub, session_id))
                        print >> f, "},"
                else:
                    session_id = func_session_path
                    print_to_file(sub,session_id)
                    fetch_path(index, os.path.join(sub, anat_session_path), os.path.join(sub, func_session_path))
                    print >> f, "},"
                    
            else:
                print "No sessions"
                session_id = ''
                print_to_file(sub,session_id)
                fetch_path(index, sub, sub)
                print >> f, "},"
        except Exception:
            print "Please make sessions are consistent across all subjects"
            raise
            
    
    try:
        print >>f, "subjects_list = ["
        for i in range(len(anat_base)):
            for sub in os.listdir(anat_base[i]): 
                #check if subject is present in subject_list
                if subject_list: 
                    if sub in subject_list:
                        print "extracting data for subject: ", sub
                        walk(i, sub)
                #check that subject is not in exclusion list
                elif sub not in exclusion_list:
                    print "extracting data for subject: ",sub
                    walk(i, sub)
    
        print >> f, "]"
        
        name = os.path.join(os.getcwd(), 'CPAC_subject_list.py')
        print "Extraction Complete...Input Subjects_list for CPAC - %s"%name
    except Exception:
        raise
    finally:
        f.close()

def generate_suplimentary_files():
    """
    Method to generate phenotypic template file
    and subject list for group analysis
    """
    from sets import Set
    import csv
    import os
    
    c = __import__('CPAC_subject_list')
    
    subject_scan_set =Set()
    subject_set = Set()
    scan_set = Set()
    data_list =[]

    for sub in c.subjects_list :
        subject_id = sub['Subject_id'] + "_" + sub['Unique_id']
        for scan in sub['rest'].keys():
            subject_scan_set.add((subject_id, scan))
            subject_set.add(subject_id)
            scan_set.add(scan)
    
    for item in subject_scan_set:
        list1 =[]
        list1.append(item[0]+"/"+item[1])
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
    if len(scan_set) >1:
        list1=['Subject_id/Scan']
        list1.extend(list(subject_set))
        list1.extend(list(scan_set))
    
    file_name = os.path.join(os.getcwd(),'phenotypic_template.csv')
    f= open(file_name, 'wb')
    writer =csv.writer(f)

    if len (scan_set) > 1:
        writer.writerow(list1)
        writer.writerows(data_list)
    else:
        writer.writerow(['Subject_id'])
        for sub in subject_set:
            writer.writerow([sub])
            
    f.close()
    
    print "Template Phenotypic file for group analysis - %s"%file_name
    
    file_name = os.path.join(os.getcwd(), "subject_list_group_analysis.txt")
    f = open(file_name, 'w')
   
    for sub in subject_set:
        print >> f, sub
    
    print "Subject list required later for group analysis - %s"%file_name
    f.close()
    
    
def run(data_config):
    """
    Run method takes data_config file as the input argument
    """
    
    import os
    import sys
    
    path, fname = os.path.split(os.path.realpath(data_config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])
    extract_data(c)
    generate_suplimentary_files()
