def merge(output_dir, scan_name, threshold, motion_f, power_f, flag):
    """
    Method to merge power parameters and motion
    parameters file
    """
    import os
    import re

    filename = scan_name + threshold + "_all_params.csv"
    filename = filename.lstrip("_")
    outfile = os.path.join(output_dir, filename)

    threshold_val = float(re.sub(r"[a-zA-Z_]", '', threshold))

    if flag:
        f = open(outfile, 'w')

        print outfile
        print >>f, "Subject,Scan,Mean_Relative_RMS_Displacement," \
        "Max_Relative_RMS_Displacement,Movements_gt_threshold,"\
        "Mean_Relative_Mean_Rotation,Mean_Relative_Maxdisp,Max_Relative_Maxdisp," \
        "Max_Abs_Maxdisp,Max_Relative_Roll,Max_Relative_Pitch," \
        "Max_Relative_Yaw,Max_Relative_dS-I,Max_Relative_dL-R," \
        "Max_Relative_dP-A,Mean_Relative_Roll,Mean_Relative_Pitch,Mean_Relative_Yaw," \
        "Mean_Relative_dS-I,Mean_Relative_dL-R,Mean_Relative_dP-A,Max_Abs_Roll," \
        "Max_Abs_Pitch,Max_Abs_Yaw,Max_Abs_dS-I,Max_Abs_dL-R,Max_Abs_dP-A," \
        "Mean_Abs_Roll,Mean_Abs_Pitch,Mean_Abs_Yaw,Mean_Abs_dS-I,Mean_Abs_dL-R,Mean_Abs_dP-A,"\
        "MeanFD,NumFD_greater_than_%.2f,rootMeanSquareFD,FDquartile(top1/4thFD),"\
        "PercentFD_greater_than_%.2f,MeanDVARS" % (threshold_val, threshold_val)
    else:
        f = open(outfile, 'a')

    motion = open(motion_f, 'r').readlines()
    power = open(power_f, 'r').readlines()

    m = motion[1].rstrip(",").split(",")
    p = power[1].rstrip(",").split(",")

    for item in m:
        f.write("%s," % item)

    for item in p[2:]:
        f.write("%s," % item)
    f.write("\n")
    f.close()


def grab(output_dir):
    """
    Method to grab all the motion parameters
    and power parameters file from each subject
    for each pipeline and merge them

    Parameters
    ----------
    output_dir : string
        Path to the datasink output directory of CPAC
    """
    import glob
    import os
    import re
    from sets import Set

    pipelines = glob.glob(os.path.join(output_dir, 'pipeline*'))

    print "number of pipelines ", len(pipelines)

    for p in pipelines:
        print "inside pipeline ->", os.path.basename(p)
        scan_list = []
        threshold_list = []

        pattern1 = re.compile(r'(\w)*scan(\w)*(\d)*(\w)*[/]')
        pattern2 = re.compile(r'(\w)*threshold_[-+]?([0-9]*\.[0-9]+|[0-9]+)')

        scans = glob.glob(os.path.join(p, '*/power_params/*/*'))

        #get the unique scans and threshold value
        for s in scans:
            val = re.search(pattern1, s)
            if val:
                scan_list.append(val.group(0).rstrip("/"))

            val = re.search(pattern2, s)
            if val:
                threshold_list.append(val.group(0))

        scan_list = Set(scan_list)
        threshold_list = Set(threshold_list)

        #print "scan_list ->", scan_list
        #print "threshold_list ->", threshold_list

        for scan in scan_list:
            for threshold in threshold_list:
                print "running for...", scan, threshold
                Flag = 1
                #merge files for each subject
                for sub  in os.listdir(p):
                    sub = os.path.join(p, sub)
                    motion_file = os.path.join(sub, 'motion_params', scan,
                                               'motion_parameters.txt')
                    power_file = os.path.join(sub, 'power_params', scan,
                                              threshold, 'pow_params.txt')

                    if os.path.exists(motion_file) and \
                        os.path.exists(power_file):
                        merge(p, scan, threshold,
                              motion_file, power_file, Flag)
                        Flag = 0

    print "process finished..."

if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        grab(sys.argv[1])
    else:
        print 'Usage: ./extract_parameters.py /path/to/datasink_dir'
