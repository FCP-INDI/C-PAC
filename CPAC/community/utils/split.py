import glob
import shutil

files = glob.glob('*.1D')
dest_list = ['group1', 'group2']
all_count = 0
bad_count = 0
not_found = []

with open('Yalegroups.txt', 'r') as grouplist:
	for line in grouplist.read().splitlines():
		if line.strip():
			if 'Group 1' in line:
				dest_flag = 0
				continue
			if 'Group 2' in line:
				dest_flag = 1
				continue
			try:
				all_count += 1
				dest_idx = files.index(line + '_rois_cc400.1D')
				print 'moving ' + files[dest_idx] + ' to ' + dest_list[dest_flag]
				shutil.move(files[dest_idx], dest_list[dest_flag])
			except ValueError:
				print line + ' not found'
				not_found.append(line)
				bad_count += 1
print '\n' + str(bad_count) + ' out of ' + str(all_count) + ' files could not be moved as they don\'t exist: \n' + str(not_found)
			
			