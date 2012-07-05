import os
import sys
import re
import commands
#import e_afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def generateMotionParameters(subject_id, scan_id, rest, movement_parameters, max_displacement):


    import os
    import numpy as np
    import re

    out_file = os.path.join(os.getcwd(), 'motion_parameters.txt')

    f= open(out_file,'w')
    print >>f, "Subject,Scan,Mean Relative RMS Displacement,"\
    "Max Relative RSM Displacement,Movements >0.1mm,Mean " \
    "Relative Mean Rotation,Mean Relative Maxdisp,Max Relative Maxdisp," \
    "Max Abs Maxdisp,Max Relative Roll,Max Relative Pitch," \
    "Max Relative Yaw,Max Relative dS-I,Max Relative dL-R," \
    "Max Relative dP-A,Mean Relative Roll,Mean Relative Pitch,Mean Relative Yaw," \
    "Mean Relative dS-I,Mean Relative dL-R,Mean Relative dP-A,Max Abs Roll," \
    "Max Abs Pitch,Max Abs Yaw,Max Abs dS-I,Max Abs dL-R,Max Abs dP-A"


    dir_name = os.path.dirname(rest)
    rest_name = os.path.basename(dir_name)
#    subject_id = os.path.basename(
#                    os.path.dirname(dir_name))
    f.write("%s," %(subject_id))
    f.write("%s," %(scan_id))
#    f.write("%s," %(rest_name))

    arr = np.genfromtxt(movement_parameters)
    arr = arr.T

    ##Relative RMS of translation
    rms= np.sqrt(arr[3]*arr[3] + arr[4]*arr[4] + arr[5]*arr[5])
    diff = np.diff(rms)
    MEANrms = np.mean(abs(diff))
    f.write("%.3f," %(MEANrms))

    MAXrms= np.max(abs(diff))
    f.write("%.3f," %(MAXrms))

    ##NUMBER OF relative RMS movements >0.1mm
    NUMmove= np.sum(abs(diff)>0.1)
    f.write("%.3f," %(NUMmove))

    ##Mean of mean relative rotation (params 1-3)
    MEANrot= np.mean(np.abs(np.diff( (abs(arr[0])+ abs(arr[1])+ abs(arr[2]))/3 ) ) )
    f.write("%.3f," %(MEANrot))

    file= open(max_displacement, 'r')
    lines =file.readlines()
    file.close()
    list1=[]

    for l in lines:
        if re.match("^\d+?\.\d+?$", l.strip()):
            list1.append(float(l.strip()))

    arr2=np.array(list1, dtype='float')

    mean=np.mean(np.diff(arr2))
    f.write("%.3f," %(mean))


    relMAX=np.max(abs(np.diff(arr2)))
    f.write("%.3f," %(relMAX))

    MAX= np.max(arr2)
    f.write("%.3f," %(MAX))


    for i in range(6):
        f.write("%.6f," %(np.max(abs(np.diff(arr[i])))))

    for i in range(6):
        f.write("%.6f," %(np.mean(np.diff(arr[i]))))

    for i in range(6):
        f.write("%.6f," %(np.max(abs(arr[i]))))

    f.close()
    return out_file


def generatePowerParams(subject_id, scan_id, rest, FD_1D, threshold, ftof_percent, sqrt_mean_raw):

    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')

    f= open(out_file,'w')
    print >>f, "Subject,Scan, MeanFD, rootMeanSquareFD, NumFD >=threshold," \
    "rmsFD, FDquartile(top 1/4th FD), PercentFD( >threshold), Num5, Num10, MeanDVARS, MeanDVARS_POW"


 #   dir_name = os.path.dirname(rest)
 #   rest_name = os.path.basename(dir_name)
 #   subject_id = os.path.basename(
 #                   os.path.dirname(dir_name))


    f.write("%s," % subject_id)
    f.write("%s," % scan_id)

    data= loadtxt(FD_1D)
    meanFD  = np.mean(data)
    f.write('%.4f,' % meanFD)

    numFD = float(data[data >=threshold].size)
    f.write('%.4f,' % numFD)

    rmsFD = np.sqrt(np.mean(data))
    f.write('%.4f,' % rmsFD)

    #Mean of the top quartile of FD is $FDquartile
    quat=int(len(data)/4)
    FDquartile=np.mean(np.sort(data)[::-1][:quat])
    f.write('%.4f,' % rmsFD)

    ##NUMBER OF FRAMES >threshold FD as percentage of total num frames
    count = np.float(data[data>threshold].size)
    percentFD = (count*100/(len(data)+1))
    f.write('%.4f,' %percentFD)


    data=loadtxt(ftof_percent)
    ###NUMBER OF relative FRAMES >5%
    num5= np.sum(data>=5)*100/len(data)
    f.write('%.4f,' % num5)

    num10= np.sum(data>=10)*100/len(data)
    f.write('%.4f,' % num10)

    meanDVARS  = np.mean(data)
    f.write('%.4f,' % meanDVARS)

    data = loadtxt(sqrt_mean_raw)
    meanDVARS_POW = np.mean(data)
    f.write('%.4f,' % meanDVARS_POW)

    f.close()
    return out_file



def get_vectors(strat):

    paths = []
    def dfs(val_list, path):

        if val_list == []:
            paths.append(list(path.split(',')[1:]))

        else:
            vals = val_list.pop()
            for val in vals:
                dfs(list(val_list), path + ',' + val)




    key_list = []
    val_list = []
    for k, v in sorted(strat.items()):
        key_list.append(k)
        val_list.append(v)

    dfs(val_list, '')

    return paths


def get_strategies(K, iterables):

    p_s_m = {}
    sm = {}
    strat = {}
    for key in K:

        vals = key.split('/')

        vals.remove('')
        for val in vals:

            if not val in p_s_m:
                p_s_m[val] = 1

    k = sorted(p_s_m.keys())
    for val in k:

        match = 0
        for variable in iterables:

            if variable in val and val.split(variable)[0] == '':

                match = 1

                it = val.split(variable)[1]
                if variable in strat:
                    item = strat[variable]
                    item.append(val)
                else:
                    strat[variable] = [val]

        if match == 0:
            strat[val] = [val]


    strat_vect = get_vectors(strat)
    return strat_vect


def get_match(sink_dir, p, strategies):

    p1 = p.split(sink_dir)[1:]
    p1 = p1[0]
    p_vals = p1.split('/')[2:]

    flag = 0
    strats = []
    for strat in strategies:
        flag = 0
        for val in p_vals:

            if not val in strat:
                flag = 1
                break
        if flag == 0:
            strats.append(str(strat))

    if strats == []:
        print "WARNING : No strategy found for ", p
        return []

    return strats

def smash(sink_dir, maps, strategies):

    new_maps = {}

    for p in sorted(maps.keys()):

        strats = get_match(sink_dir, p, strategies)

        for strat in strats:
            if strat in new_maps:
                list1 = list(new_maps[strat])
                new_maps[strat] = sorted(list(set(maps[p] + list1)))
            else:
                new_maps[strat] = sorted(list(maps[p]))

    return dict(new_maps)



def make_links(new_maps, sink_dir, sub, labels):


    print 'creating symbolink links for %s ...' % (sub)
    sym_path = os.path.join(sink_dir, 'sym_links')

    strats = sorted(new_maps.keys())
    orig_sub_dir = os.path.join(sink_dir, sub)
    wfs = sorted(os.listdir(orig_sub_dir))
    sessions = os.listdir(os.path.join(orig_sub_dir, wfs[1]))

    for strat in strats:
        l_path = os.path.join(sym_path, labels[strat])

        cmd = 'mkdir -p %s' % (l_path)
        commands.getoutput(cmd)

        sub_path = os.path.join(l_path, sub)
        cmd = 'mkdir -p %s' % (sub_path)
        commands.getoutput(cmd)


        for wf in wfs:

            for file in new_maps[strat]:

                if '/' + wf + '/' in file:
                    splits = file.split('/session_id_')
                    session = ''
                    if len(splits) == 2:
                        prepath, sessionpath = file.split('/session_id_')
                        session, postpath = (sessionpath.lstrip()).split('/', 1)
                    else:
                        session = splits[1]


                    sessionpath = os.path.join(sub_path, 'session_' + session)
                    cmd = 'mkdir -p %s' % (sessionpath)
                    commands.getoutput(cmd)

                    wf_path = os.path.join(sessionpath, wf)
                    d_path = os.path.join(sessionpath, 'Derivatives')

                    if wf in ['sca', 'vmhc', 'reho', 'alff']:
                        wf_path = d_path
                        cmd = 'mkdir -p %s' % d_path
                        wf_path = os.path.join(d_path, wf)

                    cmd = 'mkdir -p %s' % (wf_path)
                    commands.getoutput(cmd)


                    if 'seeds' in file:
                        prepath, seedpath = file.split('/seeds_')
                        seed, postpath = (seedpath.lstrip()).split('/', 1)
                        wf_path = os.path.join(wf_path, seed)
                        cmd = 'mkdir -p %s' % (wf_path)
                        commands.getoutput(cmd)

                    cmd = 'ln -s %s %s' % (file, os.path.join(wf_path, os.path.basename(file)))
                    commands.getoutput(cmd)



def make_sym_links(sink_dir, strategies, subj_list, labels):


    for subject in subj_list:
        maps = {}

        dfs_string = commands.getoutput("find %s/ -type f | perl -lne 'print tr:/::, \" $_\"' | sort -n | cut -d' ' -f2" % (os.path.join(sink_dir, subject)))

        dfs_files = dfs_string.split('\n')

        for file in dfs_files:
            pdir = os.path.dirname(file)

            val = None
            if pdir in maps:
                val = maps[pdir]
                val.append(file)
            else:
                val = [file]
            maps[os.path.dirname(file)] = val

        new_maps = smash(sink_dir, maps, strategies)

        make_links(new_maps, sink_dir, subject, labels)


def write_label_linkage(sink_dir, strategies, proper_names):

    labels = {}

    sym_path = os.path.join(sink_dir, 'sym_links')
    cmd = 'mkdir -p %s' % (sym_path)
    commands.getoutput(cmd)

    f = open('%s/label_linkage.txt' % (sym_path), 'w')
    label_idx = 1
    for strategy in strategies:

        labels[str(strategy)] = 'Strategy_' + str(label_idx)
        label_idx += 1
        print >>f, '\n\n', labels[str(strategy)]
        print >>f, '___________\n'
        for parameter in strategy:

            (strat, value) = parameter.rsplit('_', 1)

            try:
                print >>f, proper_names[strat], value
            except:
                pass

    f.close()
    return labels


def symlink_creator(sink_dir):

    import os
    subj_list = sorted(os.listdir(sink_dir))
    subjectDir = os.path.join(sink_dir, subj_list[0])

    iterables = ['threshold', 'csf_threshold', 'fwhm', 'gm_threshold', 'hp', 'lp', 'nc', 'scrubbed', 'selector', 'target_angle', 'wm_threshold']
    proper_names = {'threshold':'Scrubbing Threshold = ', 'csf_threshold':'Cerebral Spinal Fluid Threshold = ', 'fwhm':'Full Width Half Max = ',
                    'gm_threshold':'Gray Matter Threshold = ', 'hp':'Temporal Filtering HighPass LowCutoff = ', 'lp':'Temporal Filtering LowPass HighCutoff = ',
                    'nc':'Compcor: Number Of Components = ', 'scrubbed':'Data Scrubbed = ', 'selector':'Nuisance Signal Corrections = ',
                    'target_angle_deg':'Median Angle Correction: Traget Angle in Degree = ', 'wm_threshold':'White Matter Threshold = '}

    #iterables = ['threshold', 'csf_threshold', 'fwhm', 'gm_threshold', 'hp', 'lp', 'nc', 'scrubbed', 'seeds', 'selector', 'session_id', 'target_angle', 'wm_threshold']
    wfs = os.listdir(subjectDir)

    labels = {}
    label_id = 0

    dfs_string = commands.getoutput("find %s/ -type f | perl -lne 'print tr:/::, \" $_\"' | sort -n | cut -d' ' -f2" % (subjectDir))

    dfs_files = dfs_string.split('\n')

    maps = {}

    m1 = {}
    for file in dfs_files:
        f = str(file)
        spath = os.path.join(sink_dir, subj_list[0])
        vmhc_path = os.path.join(spath, 'vmhc')
        alff_path = os.path.join(spath, 'alff')
        nuisance_path = os.path.join(spath, 'nuisance')
        sca_path = os.path.join(spath, 'sca')
        scrubbing_path = os.path.join(spath, 'scrubbing')
        anat_path = os.path.join(spath, 'anat')
        func_path = os.path.join(spath, 'func')
        seg_path = os.path.join(spath, 'segment')
        reg_path = os.path.join(spath, 'reg')
        param_path = os.path.join(spath, 'parameters')

        file = re.sub(r'%s' % (vmhc_path), '', file)
        file = re.sub(r'%s' % (alff_path), '', file)
        file = re.sub(r'%s' % (nuisance_path), '', file)
        file = re.sub(r'%s' % (sca_path), '', file)
        file = re.sub(r'%s' % (scrubbing_path), '', file)
        file = re.sub(r'%s' % (anat_path), '', file)
        file = re.sub(r'%s' % (func_path), '', file)
        file = re.sub(r'%s' % (reg_path), '', file)
        file = re.sub(r'%s' % (seg_path), '', file)
        file = re.sub(r'%s' % (param_path), '', file)

        #print file
        pdir = os.path.dirname(file)

        val = None
        v = None
        if pdir in maps:
            val = maps[pdir]
            v = m1[pdir]
            val.append(file)
            v.append(f)
        else:
            val = [file]
            v = [f]
        maps[os.path.dirname(file)] = val
        m1[os.path.dirname(file)] = v

    strategies = sorted(list(get_strategies(maps.keys(), iterables)))

    labels = write_label_linkage(sink_dir, strategies, proper_names)

    make_sym_links(sink_dir, strategies, subj_list, labels)


def getIndx(in_file):

       indx = []

       if(isinstance(in_file, list)):
               for file in in_file:
                       f = open(file, 'r')
                       line = f.readline()
                       line = line.strip(',')
                       indx.append(map(int, line.split(",")))
                       f.close()
               print "indx ", indx
               return indx
       else:
            f = open(file, 'r')
            line = f.readline()
            line = line.strip(',')
            indx.append(map(int, line.split(",")))
            f.close()
            print "indx in else", indx
            return indx


def getStartIdx(in_file):

    start_indx = []

    if(isinstance(in_file, list)):
        for file in in_file:
            f = open(file, 'r')
            line = f.readline()
            start_indx.append(int(line.split(",")[0]))
            f.close()
        return start_indx
    else:
        f = open(in_file, 'r')
        myList = []
        line = f.readline()
        myList.extend(line.split(","))
        f.close()
        return int(myList[0])


def getStopIdx(in_file):

    stop_indx = []
    if (isinstance(in_file, list)):
        for file in in_file:
            f = open(file, 'r')
            line = f.readline()
            myList = []
            myList.extend(line.split(","))
            length = len(myList)
            f.close()

        if str(myList[length - 1]) == "":
            stop_indx.append(
                int(myList[length - 2]))
        else:
            stop_indx.append(
                int(myList[length - 1]))
        return stop_indx
    else:
        f = open(in_file, 'r')
        line = f.readline()
        myList = []
        myList.extend(line.split(","))
        length = len(myList)
        f.close()
        if str(myList[length - 1]) == "":
            return int(myList[length - 2])
        else:
            return int(myList[length - 1])


def last_vol(nvols, stopIdx, startIdx):

    vol = nvols
    last_volume = (int(vol) - 1)
    return last_volume


def TRendminus1(nvols, stopIdx, startIdx):

    vol = nvols
    last_vol_minus_one = None

    last_volume = (int(vol) - 1)
    last_vol_minus_one = last_volume - 1

    return last_vol_minus_one



def scCopy(in_file):

    import os

    cwd = os.getcwd()

    out_file = 'pow_params.txt'
    path = os.path.join(cwd, out_file)
    out_file = path

    dir_name = os.path.dirname(in_file)
    rest_name = os.path.basename(dir_name)
    subject_id = os.path.basename(
                    os.path.dirname(dir_name))

    f = open(out_file, 'a')

    f.write("%s," % subject_id)
    f.write("%s," % rest_name)

    f.close()

    return out_file


def createSC(in_file):
    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'FD.1D')

    cmd1 = sb.Popen(
        ['awk', '{x=$4} {y=$5} {z=$6} {a=$1} {b=$2} {c=$3} ' +
        '{print 2*3.142*50*(a/360),2*3.142*50*(b/360),' +
        ' 2*3.142*50*(c/360), x, y, z}',
        in_file], stdin=sb.PIPE, stdout=sb.PIPE,)

    cmd2 = sb.Popen(
        ['awk', '{a=$1} {b=$2} {c=$3} {x=$4} {y=$5} {z=$6} ' +
        'NR>=1{print a-d, b-e, c-f, x-u, y-v, z-w}' +
        '{d=a} {e=b} {f=c} {u=x} {v=y} {w=z}'],
        stdin=cmd1.stdout, stdout=sb.PIPE,)
    cmd3 = sb.Popen(
        ['awk', '{ for (i=1; i<=NF; i=i+1) {' +
        'if ($i < 0) $i = -$i} print}'],
        stdin=cmd2.stdout, stdout=sb.PIPE,)
    cmd4 = sb.Popen(
        ['awk', '{a=$1+$2+$3+$4+$5+$6} {print a}'],
        stdin=cmd3.stdout, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd4.communicate()
    output = stdout_value.split()
    f = open(out_file, 'w')
    for out in output:
        print >>f, float(out)
    f.close()
    return out_file


def setMeanFD(infile_a, infile_b):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    data = loadtxt(infile_b)
    mean = np.mean(data)
    print "mean ->", mean

    f = open(out_file, 'a')
    f.write('%.4f,' % mean)

    f.close()

    return out_file


def setNumFD(infile_a, infile_b, threshold):

    import subprocess as sb
    import os
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    data = loadtxt(infile_b)
    count = float(data[data >= threshold].size)
    print "count ->", count

    f = open(out_file, 'a')
    f.write('%.4f,' % count)

    f.close()

    return out_file


def setPercentFD(infile_a, infile_b, threshold):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    data = loadtxt(infile_b)
    count = np.float(data[data > threshold].size)
    sum = (count*100/(len(data)+1))
    print "sum ->", sum

    f = open(out_file, 'a')
    f.write('%.4f,' %sum)

    f.close()

    return out_file


def setFramesEx(in_file, threshold):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'frames_ex.1D')
    data = loadtxt(in_file)
    #masking zeroth timepoint value as 0, since the mean displacment value for
    #zeroth timepoint cannot be calculated, as there is no timepoint before it
    data[0] = 0

    extra_indices = []

    indices = [i[0] for i in (np.argwhere(data >= threshold)).tolist()]

    #adding addtional 2 after and one before timepoints
    for i in indices:

     if i > 0:
        extra_indices.append(i-1)
        if i+1 < data.size:
            extra_indices.append(i+1)
        if i+2 < data.size:
            extra_indices.append(i+2)

    indices = list(set(indices) | set(extra_indices))


    f = open(out_file, 'a')

    print "indices->", indices
    for idx in indices:
        f.write('%s,' % int(idx))

    f.close()

    return out_file

def setFramesIN(in_file, threshold, exclude_list):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'frames_in.1D')

    data = loadtxt(in_file)
    #masking zeroth timepoint value as 0, since the mean displacment value for
    #zeroth timepoint cannot be calculated, as there is no timepoint before it
    data[0] = 0

    indices = [i[0] for i in (np.argwhere(data < threshold)).tolist()]

    indx = []
    f = open(exclude_list, 'r')
    line = f.readline()
    if line:
        line = line.strip(',')
        indx = map(int, line.split(","))
    f.close()
    print indx

    if indx:
        indices = list(set(indices) - set(indx))

    f = open(out_file, 'a')

    for idx in indices:
        f.write('%s,' % int(idx))

    f.close()

    return out_file

def setFramesInList(in_file):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'frames_in_list.1D')

    indices = []
    f = open(in_file, 'r')
    line = f.readline()
    if line:
        line = line.strip(',')
        indices = map(int, line.split(","))
    f.close()

    f = open(out_file, 'a')

    for idx in indices:
        print >>f, int(idx)

    f.close()

    return out_file


def setMeanDVARS(infile_a, infile_b):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    data = loadtxt(infile_b)
    mean = np.mean(data)

    f = open(out_file, 'a')

    f.write('%.4f,' % mean)

    f.close()

    return out_file


def setNUM5(infile_a, infile_b):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    ###NUMBER OF relative FRAMES >5%
    frame_percentage = 5

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    data = loadtxt(infile_b)
    count = np.float(data[data >= frame_percentage].size)

    f = open(out_file, 'a')

    f.write('%.4f,' % count)

    f.close()

    return out_file


def setNUM10(infile_a, infile_b):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    ###NUMBER OF relative FRAMES >10%
    frame_percentage = 10

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    copycmd = sb.Popen(['cp', infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE )
    out_val, error_val = copycmd.communicate()
    print "outval ---> ", out_val
    print "error_val -->", error_val

    data = loadtxt(infile_b)
    count = np.float(data[data >= frame_percentage].size)

    f = open(out_file, 'a')
    f.write('%.4f,' % count)

    f.close()

    return out_file


def setSqrtMeanDeriv(in_file):
    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'sqrt_mean_deriv_sq.1D')

    cmd = sb.Popen(['awk', '{x=$1} {printf( "%.6f ", sqrt(x))}',
                    in_file], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
        print >>f, float(out)

    f.close()

    return out_file


def setSqrtMeanRaw(in_file):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'sqrt_mean_raw_sq.1D')

    cmd = sb.Popen(['awk', '{x=$1} {printf( "%.6f ", sqrt(x))}',
                    in_file], stdin=sb.PIPE, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
        print >>f, float(out)

    f.close()

    return out_file


def setFtoFPercentChange(infile_a, infile_b):

    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'ftof_percent_change.1D')

    cmd1 = sb.Popen(['awk', 'NR==FNR{a[NR]=$1; next} {print a[FNR], $1}',
                     infile_a, infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    cmd2 = sb.Popen(['awk', '{x=$1} {y=$2} {printf( "%.6f ", ((x/y)*100))}'],
                    stdin=cmd1.stdout, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd2.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
        print >>f, float(out)

    f.close()

    return out_file


def setNUMFD(in_file, threshold):

    import subprocess as sb
    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'numFD')

    data = loadtxt(in_file)
    count = np.float(data[data >= threshold].size)
    print "count ->", count

    f = open(out_file, 'a')
    f.write('%.4f,' % count)

    f.close()

    return out_file



def setScrubbedMotion(infile_a, infile_b):

    import os

    out_file = os.path.join(os.getcwd(), 'rest_mc_scrubbed.1D')

    f1= open(infile_a)
    f2=open(infile_b)
    l1=f1.readline()
    l2=f2.readlines()
    f1.close()
    f2.close()


    l1=l1.rstrip(',').split(',')

    f = open(out_file, 'a')
    for l in l1:
        data=l2[int(l.strip())]
        f.write(data)
    f.close()
    return out_file


def pToFile(time_series):

    import os
    import re
    import commands
    import sys

    dir = os.path.dirname(time_series)




    dir1 = re.sub(r'(.)+_seeds_', '', dir)

    dir1 = dir1.split('.nii.gz')
    dir1 = dir1[0]

    ts_oneD = os.path.join(dir, dir1 + '.1D')
    cmd = "cp %s %s" % (time_series, ts_oneD)
    print cmd

    sys.stderr.write('\n' + commands.getoutput(cmd))
    return os.path.abspath(ts_oneD)


def mean_roi_signal(data_volume, roi_mask):
    import numpy as np
    Y = data_volume[roi_mask].T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    return Yc.mean(1)


def get_standard_background_img(in_file, file_parameters):
    import os
    from nibabel import load
    img = load(in_file)
    hdr = img.get_header()
    group_mm = int(hdr.get_zooms()[2])
    FSLDIR, MNI = file_parameters
    print "gorup_mm -> ", group_mm
    path = FSLDIR + '/data/standard/' + MNI +'_T1_%smm_brain.nii.gz' % (group_mm)
    print "path ->", path
    return os.path.abspath(path)

def get_nvols(in_file):

    from nibabel import load
    img = load(in_file)
    hdr = img.get_header()
    n_vol = int(hdr.get_data_shape()[3])
    op_string = '-abs -bin -Tmean -mul %d' % (n_vol)
    return op_string


def copyGeom(infile_a, infile_b):
    import subprocess as sb
    out_file = infile_b
    cmd = sb.Popen(['fslcpgeom',
                   infile_a, out_file], stdin=sb.PIPE, stdout=sb.PIPE,)
    stdout_value, stderr_value = cmd.communicate()
    return out_file


def getTuple(infile_a, infile_b):

    return (infile_a, infile_b[1])


def getOpString(mean, std_dev):

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))

    op_string = str1 + " -mas %s"

    return op_string


def pick_wm_0(probability_maps):

    import sys
    import os

    if(isinstance(probability_maps, list)):

        if(len(probability_maps) == 1):
            probability_maps = probability_maps[0]
        for file in probability_maps:
            print file
            if file.endswith("prob_0.nii.gz"):

                return file
    return None


def pick_wm_1(probability_maps):

    import sys
    import os

    if(isinstance(probability_maps, list)):

        if(len(probability_maps) == 1):
            probability_maps = probability_maps[0]
        for file in probability_maps:
            print file
            if file.endswith("prob_1.nii.gz"):

                return file
    return None


def pick_wm_2(probability_maps):

    import sys
    import os
    if(isinstance(probability_maps, list)):

        if(len(probability_maps) == 1):
            probability_maps = probability_maps[0]
        for file in probability_maps:
            print file
            if file.endswith("prob_2.nii.gz"):

                return file
    return None

def get_idx(in_files, stop_idx, start_idx):

    stopidx = None
    startidx = None
    from nibabel import load

    img = load(in_files)
    hdr = img.get_header()
    nvols = int(hdr.get_data_shape()[3])

    if (start_idx == None) or (start_idx < 0) or (start_idx > (nvols - 1)):
        startidx = 0
    else:
        startidx = start_idx

    if (stop_idx == None) or (stop_idx > (nvols - 1)):
        stopidx = nvols - 1
    else:
        stopidx = stop_idx

    return stopidx, startidx


def getImgNVols(in_files, stopIdx, startIdx):

    out = []
    from nibabel import load
    if(isinstance(in_files, list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()

            if len(hdr.get_data_shape()) > 3:
                nvols = int(hdr.get_data_shape()[3])
            else:
                nvols = 1
            out.append(nvols)
        return out

    else:
        img = load(in_files)
        hdr = img.get_header()
        if len(hdr.get_data_shape()) > 3:
            nvols = int(hdr.get_data_shape()[3])
        else:
            nvols = 1
        return [nvols]



def getImgTR(in_files, TRa):

    out = []
    from nibabel import load
    if(isinstance(in_files, list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()
            tr = float(hdr.get_zooms()[3])
            if tr > 10:
                tr = float(float(tr) / 1000.0)
            if not (TRa == None):
                diff = None
                if TRa > tr:
                    diff = TRa - tr
                else:
                    diff = tr - TRa

                if (diff > 0.001):
                    print "Warning: specified TR  %f and TR in image header  %f do not match:" % (TRa, tr)
                out.append(TRa)
            else:
                out.append(tr)
        return out
    else:
        img = load(in_files)
        hdr = img.get_header()
        tr = float(hdr.get_zooms()[3])
        if tr > 10:
            tr = float(float(tr) / 1000.0)
        if not (TRa == None):
            diff = None
            if TRa > tr:
                diff = TRa - tr
            else:
                diff = tr - TRa

            if (diff > 0.001):
                print "Warning: specified TR  %f and TR in image header  %f do not match:" % (TRa, tr)
            return [TRa]
        else:
            return [tr]


def getN1(TR, nvols, HP):

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n1 = int("%1.0f" % (float(n_lp - 1.0)))

    return n1


def getN2(TR, nvols, LP, HP):

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n_hp = float(LP) * float(int(nvols)) * float(TR)
    n2 = int("%1.0f" % (float(n_hp - n_lp + 1.0)))

    return n2


def takemod(nvols):

    decisions = []
    for vol in nvols:
        mod = int(int(vol) % 2)

        if mod == 1:
            decisions.append([0])
            #return [0]
        else:
            decisions.append([1])
            #return [1]

    return decisions


def set_op_str(n2):

    strs = []
    for n in n2:
        str = "-Tmean -mul %f" % (n)
        strs.append(str)
    return strs


def set_op1_str(nvols):

    strs = []
    for vol in nvols:
        str = '-Tmean -mul %d -div 2' % (int(vol))
        strs.append(str)

    return strs


def set_gauss(fwhm):

    op_string = ""

    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string


def getEXP(nvols):

    expr = []
    for vol in nvols:
        vol = int(vol)
        expr.append("'a*sqrt('%d'-3)'" % vol)

    return expr


def select(pp, scrubbed, sc_pp):

    pp_s = None
    if (scrubbed):

        pp_s = sc_pp

    else:

        pp_s = pp

    return pp_s


def selectM(movement_parameters,
            scrubbed, scrubbed_movement_parameters):

    mp = None
    if (scrubbed):

        mp = scrubbed_movement_parameters

    else:

        mp = movement_parameters

    return mp



def selector_wf():

    swf = pe.Workflow(name='scrubbing_selector')

    inputNode = pe.Node(util.IdentityInterface(fields=[
                                            'preprocessed',
                                            'movement_parameters',
                                            'scrubbed_preprocessed',
                                            'scrubbed_movement_parameters'
                                            ]),
                        name='inputspec')

    iRun = pe.Node(util.IdentityInterface(fields=['scrubbed']),
                             name='scrubbed_input')

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                'preprocessed_selector',
                                                'movement_parameters_selector'
                                                ]),
                        name='outputspec')

    selectP = pe.MapNode(util.Function(input_names=[
                                              'pp',
                                              'scrubbed',
                                              'sc_pp'],
                                      output_names=['pp_s'],
                        function=select),
                        name='selectP', iterfield=['pp',
                                                'sc_pp'])

    selectP1 = pe.MapNode(util.Function(input_names=[
                                        'movement_parameters',
                                        'scrubbed',
                                        'scrubbed_movement_parameters'],
                                        output_names=['mp'],
                          function=selectM),
                          name='selectP1',
                          iterfield=[
                                'movement_parameters',
                                'scrubbed_movement_parameters'])

    swf.connect(inputNode, 'preprocessed',
                selectP, 'pp')
    swf.connect(inputNode, 'scrubbed_preprocessed',
                selectP, 'sc_pp')
    swf.connect(iRun, 'scrubbed',
                selectP, 'scrubbed')

    swf.connect(inputNode, 'movement_parameters',
                selectP1, 'movement_parameters')
    swf.connect(inputNode, 'scrubbed_movement_parameters',
                selectP1, 'scrubbed_movement_parameters')
    swf.connect(iRun, 'scrubbed',
                selectP1, 'scrubbed')

    swf.connect(selectP, 'pp_s',
                outputNode, 'preprocessed_selector')
    swf.connect(selectP1, 'mp',
                outputNode, 'movement_parameters_selector')

    return swf


def getStatsDir(in_files):

    import os

    stats_dir = []

    for in_file in in_files:

        parent_path = os.path.dirname(in_file)

        stats_dir.append(parent_path + '/stats')

    return stats_dir


def create_alff_dataflow(name,
                         sublist,
                         analysisdirectory,
                         rest_name,
                         at,
                         atw):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(
                                    infields=['subject_id'],
                                    outfields=['rest_res',
                                        'rest_mask',
                                        'rest_mask2standard']),
                         name=name)
    datasource.inputs.base_directory = analysisdirectory
    datasource.inputs.template = at
    datasource.inputs.template_args = dict(
                                    rest_res=[['subject_id',
                                              rest_name + '_res']],
                                    rest_mask=[['subject_id',
                                              rest_name + '_mask']],
                                    rest_mask2standard=[['subject_id',
                                    rest_name + '_mask2standard']])
    datasource.iterables = ('subject_id', sublist)

    datasource_warp = pe.Node(interface=nio.DataGrabber(
                                            infields=['subject_id'],
                                            outfields=['premat',
                                                'fieldcoeff_file']),
                              name=name + '_warp')

    datasource_warp.inputs.base_directory = analysisdirectory
    datasource_warp.inputs.template = atw
    datasource_warp.inputs.template_args = dict(
                                    premat=[['subject_id',
                                            'example_func2highres.mat']],
                                    fieldcoeff_file=[['subject_id',
                                        'highres2standard_warp.nii.gz']])
    datasource_warp.iterables = ('subject_id',
                                sublist)

    return datasource, datasource_warp


def create_sca_dataflow(name, sublist, analysisdirectory, rt, rtw):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(
                                    infields=['subject_id'],
                                    outfields=['rest_res2standard',
                                               'rest_res_filt',
                                               'rest_mask2standard',
                                               'ref']),
                         name=name)
    datasource.inputs.base_directory = analysisdirectory
    #datasource.inputs.template = '%s/*/%s.nii.gz'
    datasource.inputs.template = rt
    datasource.inputs.template_args = dict(
                        rest_res2standard=[['subject_id',
                                rest_name + '_res2standard']],
                        rest_res_filt=[['subject_id',
                                rest_name + '_res_filt']],
                        rest_mask2standard=[['subject_id',
                                rest_name + '_mask2standard']],
                        ref=[['subject_id',
                                'example_func']])
    datasource.iterables = ('subject_id',
                            sublist)

    datasource_warp = pe.Node(interface=nio.DataGrabber(
                                        infields=['subject_id'],
                                        outfields=['premat',
                                                   'fieldcoeff_file',
                                                   'warp',
                                                   'postmat']),
                              name=name + '_warp')
    datasource_warp.inputs.base_directory = analysisdirectory
    #datasource_warp.inputs.template = '%s/*/%s.nii.gz'
    datasource_warp.inputs.template = rtw
    datasource_warp.inputs.template_args = dict(
                                        premat=[['subject_id',
                                        'example_func2highres.mat']],
                                        fieldcoeff_file=[['subject_id',
                                        'highres2standard_warp.nii.gz']],
                                        warp=[['subject_id',
                                        'stand2highres_warp.nii.gz']],
                                        postmat=[['subject_id',
                                        'highres2example_func.mat']])
    datasource_warp.iterables = ('subject_id', sublist)

    return datasource, datasource_warp


def create_seed_dataflow(seed_file):

    import nipype.interfaces.io as nio
    import os
    import re
    f = open(seed_file, 'r')

    seeds = f.readlines()
    seed_dir = None

    seed_list = []
    seed1 = (seeds[0]).rstrip('\r\n')
    seed1 = (seed1.split(os.path.basename(seed1)))[0]
    seed1 = re.sub('home/', 'home2/', seed1)
    seed_dir = seed1
    for seed in seeds:
        seed_path = seed.rstrip('\r\n')
        seed_list.append(os.path.basename(seed_path))

    f.close()

    datasource = pe.Node(interface=nio.DataGrabber(
                                        infields=['seeds'],
                                        outfields=['out_file']),
                         name="datasource_seeds")
    datasource.inputs.base_directory = seed_dir
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(out_file='%s')
    datasource.inputs.template_args = dict(out_file=[['seeds']])
    datasource.iterables = ('seeds', seed_list)

    return datasource


def create_vmhc_dataflow(name,
                         sublist,
                         analysisdirectory,
                         anat_name,
                         rest_name,
                         vt,
                         vta,
                         vtw):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    datasource = pe.Node(interface=nio.DataGrabber(
                                    infields=['subject_id'],
                                    outfields=['rest_res']),
                         name=name + '_res')
    datasource.inputs.base_directory = analysisdirectory
    datasource.inputs.template = vt
    datasource.inputs.template_args = dict(
                                    rest_res=[['subject_id',
                                    rest_name + '_res']])
    datasource.iterables = ('subject_id',
                            sublist)

    da = pe.Node(interface=nio.DataGrabber(
                                    infields=['subject_id'],
                                    outfields=['reorient',
                                               'brain']),
                 name=name + '_RPI')
    da.inputs.base_directory = analysisdirectory
    da.inputs.template = vta
    da.inputs.template_args = dict(reorient=[['subject_id',
                                    anat_name + '_RPI']],
                                   brain=[['subject_id',
                                    anat_name + '_brain']])
    da.iterables = ('subject_id',
                   sublist)

    datasource_warp = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                outfields=['example_func2highres_mat']),
                              name=name + '_example_func')
    datasource_warp.inputs.base_directory = analysisdirectory
    #datasource_warp.inputs.template = '%s/*/%s.nii.gz'
    datasource_warp.inputs.template = vtw
    datasource_warp.inputs.template_args = dict(example_func2highres_mat=[['subject_id',
                                                'example_func2highres.mat']])
    datasource_warp.iterables = ('subject_id',
                                sublist)

    return datasource, da, datasource_warp


def create_gp_dataflow(base_dir, modelist, dervlist, labelist, template_dict, template_args, name):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces.utility import IdentityInterface


    wf = pe.Workflow(name='gp_dataflow')


    inputnode = pe.Node(interface=IdentityInterface(fields=['model', 'derivative', 'label'],
                                                    mandatory_inputs=True),
                        name='inputnode')
    inputnode.iterables = [('model', modelist),
                           ('derivative', dervlist),
                           ('label', labelist)]

    datasource = pe.Node(interface=nio.DataGrabber(infields=['model_name', 'derivative', 'label'],
                                                   outfields=['mat', 'con', 'fts', 'grp', 'derv']),
                         name=name)
    datasource.inputs.base_directory = base_dir
    datasource.inputs.template = '*'

    datasource.inputs.field_template = template_dict
    datasource.inputs.template_args = template_args

    wf.connect(inputnode, 'model', datasource, 'model_name')
    wf.connect(inputnode, 'derivative', datasource, 'derivative')
    wf.connect(inputnode, 'label', datasource, 'label')
    return wf


def extract_subjectID(out_file):

    import sys

    outs = (out_file.split('_subject_id_'))[1]
    out_file = (outs.split('/'))[0]
#    print "extract_subjectID out_file", out_file
#    sys.stderr.write('\n>>>>>D<<<<< ' + out_file)

    return out_file


def listdir_nohidden(path):
    import os
    list1 = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            list1.append(os.path.splitext(os.path.splitext(f)[0])[0])
    return list1


def create_parc_dataflow(unitDefinitionsDirectory):

    import nipype.interfaces.io as nio
    import os

    unitlist = [os.path.splitext(os.path.splitext(f)[0])[0]
                for f in os.listdir(unitDefinitionsDirectory)]
#    print "unitList ->", unitlist
    datasource = pe.Node(interface=nio.DataGrabber(
                                infields=['parcelation'],
                                outfields=['out_file']),
                         name="datasource_parc")
    datasource.inputs.base_directory = unitDefinitionsDirectory
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(out_file='%s.nii.gz')
    datasource.inputs.template_args = dict(out_file=[['parcelation']])
    datasource.iterables = ('parcelation',
                            unitlist)

    return datasource


def create_mask_dataflow(voxelMasksDirectory):

    import nipype.interfaces.io as nio
    import os

    masklist = [os.path.splitext(
                    os.path.splitext(f)[0])[0]
                    for f in os.listdir(voxelMasksDirectory)]
#    print masklist
    datasource = pe.Node(interface=nio.DataGrabber(
                                        infields=['mask'],
                                        outfields=['out_file']),
                         name="datasource_mask")
    datasource.inputs.base_directory = voxelMasksDirectory
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(out_file='%s.nii.gz')
    datasource.inputs.template_args = dict(out_file=[['mask']])
    datasource.iterables = ('mask', masklist)

    return datasource


def gen_csv_for_parcelation(data_file,
                            template,
                            unitTSOutputs):

    import nibabel as nib
    import csv
    import numpy as np
    import os

    unit = nib.load(template)
    unit_data = unit.get_data()
    datafile = nib.load(data_file)
    img_data = datafile.get_data()
    vol = img_data.shape[3]

    nodes = np.unique(unit_data).tolist()
    sorted_list = []
    node_dict = {}
    out_list = []

    for n in nodes:
        if n > 0:
            node_array = img_data[unit_data == n]
            node_array = node_array.T
            time_points, no_of_voxels = node_array.shape
            list1 = [n]
            node_str = 'node %s' % (n)
            node_dict[node_str] = node_array
            for t in range(0, time_points):
                avg = node_array[t].mean()
                avg = np.round(avg, 3)
                list1.append(avg)
            sorted_list.append(list1)

    tmp_file = os.path.splitext(
                    os.path.basename(template))[0]
    tmp_file = os.path.splitext(tmp_file)[0]
    csv_file = os.path.abspath('parc_' + tmp_file + '.csv')
    numpy_file = os.path.abspath('parc_' + tmp_file + '.npz')

    if unitTSOutputs[0]:
        f = open(csv_file, 'wt')
        writer = csv.writer(f, delimiter=',',
                                quoting=csv.QUOTE_MINIMAL)
        headers = ['node/volume']
        for i in range(0, vol):
            headers.append(i)
        writer.writerow(headers)
        writer.writerows(sorted_list)
        f.close()
        out_list.append(csv_file)

    if unitTSOutputs[1]:
        np.savez(numpy_file, **dict(node_dict))
        out_list.append(numpy_file)

    return out_list


def gen_csv_for_mask(data_file,
                     template,
                     voxelTSOutputs):

    import nibabel as nib
    import numpy as np
    import csv
    import os

    unit = nib.load(template)
    unit_data = unit.get_data()
    datafile = nib.load(data_file)
    img_data = datafile.get_data()
    header_data = datafile.get_header()
    qform = header_data.get_qform()
    sorted_list = []
    vol_dict = {}
    out_list = []

    node_array = img_data[unit_data != 0]
    node_array = node_array.T
    time_points = node_array.shape[0]
    for t in range(0, time_points):
        str = 'vol %s' % (t)
        vol_dict[str] = node_array[t]
        val = node_array[t].tolist()
        val.insert(0, t)
        sorted_list.append(val)

    tmp_file = os.path.splitext(
                  os.path.basename(template))[0]
    tmp_file = os.path.splitext(tmp_file)[0]
    csv_file = os.path.abspath('mask_' + tmp_file + '.csv')
    numpy_file = os.path.abspath('mask_' + tmp_file + '.npz')

    if voxelTSOutputs[0]:
        f = open(csv_file, 'wt')
        writer = csv.writer(f, delimiter=',',
                            quoting=csv.QUOTE_MINIMAL)
        one = np.array([1])
        headers = ['volume/xyz']
        cordinates = np.argwhere(unit_data != 0)
        for val in range(np.alen(cordinates)):
            ijk_mat = np.concatenate([cordinates[val], one])
            ijk_mat = ijk_mat.T
            product = np.dot(qform, ijk_mat)
            val = tuple(product.tolist()[0:3])
            headers.append(val)
        writer.writerow(headers)
        writer.writerows(sorted_list)
        f.close()
        out_list.append(csv_file)

    if voxelTSOutputs[1]:
        np.savez(numpy_file, **dict(vol_dict))
        out_list.append(numpy_file)

    return out_list


def gen_csv_for_surface(rh_surface_file,
                        lh_surface_file,
                        verticesTSOutputs):

    import gradunwarp
    import numpy as np
    import os
    out_list = []

    if verticesTSOutputs[0]:
        rh_file = os.path.splitext(
                        os.path.basename(rh_surface_file))[0] + '_rh.csv'
        mghobj1 = gradunwarp.mgh.MGH()

        mghobj1.load(rh_surface_file)
        vol = mghobj1.vol
        (x, y) = vol.shape
#        print "rh shape", x, y

        np.savetxt(rh_file, vol, delimiter=',')
        out_list.append(rh_file)

        lh_file = os.path.splitext(os.path.basename(lh_surface_file))[0] + '_lh.csv'
        mghobj2 = gradunwarp.mgh.MGH()

        mghobj2.load(lh_surface_file)
        vol = mghobj2.vol
        (x, y) = vol.shape
#        print "lh shape", x, y

        np.savetxt(lh_file,
                   vol,
                   delimiter=',')
        out_list.append(lh_file)

    return out_list


def create_datasink(source_dir):

    import nipype.interfaces.io as nio
    print"create data sink"
    datasink = pe.Node(nio.DataSink(), name='data_sink')
    datasink.inputs.base_directory = source_dir
    datasink.inputs.regexp_substitutions = [(r"[/](_)+", '/'), (r"^(_)+", ''), (r"[\w]*rename(\d)+[/]", ''),
                                            (r"_subject_(\w)*(\d)*", ''), (r"apply_warp(\d)+[/]", 'MNI_Outputs/'),
                                            (r"func_(\w)+(\d)+[/]", ''), (r"median_angle(\d)+", ''),
                                            (r"sc_(\d)*(\w)+(\d)+[/]", ''), (r"seg_(\w)+(\d)+[/]", ''),
                                            (r"timeseries_(\w)+(\d)+[/]", ''), (r"(\.nii.gz)[/]", '/')]
    return datasink

def create_gp_datasink(source_dir):

    import nipype.interfaces.io as nio
    import os
    print"create gp data sink"
    datasink = pe.Node(nio.DataSink(), name='gp_data_sink')
    dir = os.path.join(source_dir, 'groupAnalysis')
    datasink.inputs.base_directory = dir
    datasink.inputs.regexp_substitutions = [(r"[/](_)+", '/'),
                                            (r"^(_)+", ''),
                                            (r"[\w]*rename(\d)+[/]", ''),
                                            (r"_label_(\w)*(\d)*(\w)*_model", '/model'),
                                            (r"gp_(\w)+(\d)+", ''),
                                            (r"Merged/model_(\w)*(\d)*(\w)*[/]", 'Merged/')]
    return datasink
