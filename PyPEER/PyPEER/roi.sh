data_dir='/Users/jakeson/Documents/CMI/peer/data'
eye_mask='/Users/jakeson/Documents/CMI/peer/MNI152_T1_2mm_eye_mask.nii.gz'

for participant in $(ls $data_dir); do

  for scan in $(ls $data_dir/$participant | grep .nii.gz); do

    echo "Running fslutils on $participant/$scan."

    fslmaths $data_dir/$participant/$scan -mas $eye_mask $data_dir/$participant/em_$scan
    fslmaths $data_dir/$participant/em_$scan -Tmean $data_dir/$participant/mean_$scan
    fslmaths $data_dir/$participant/em_$scan -Tstd $data_dir/$participant/std_$scan
    fslmaths $data_dir/$participant/em_$scan -sub $data_dir/$participant/mean_$scan -div $data_dir/$participant/std_$scan $data_dir/$participant/roi_$scan

    rm $data_dir/$participant/em_$scan
    rm $data_dir/$participant/mean_$scan
    rm $data_dir/$participant/std_$scan

    echo "$participant/$scan roi saved."

  done

done
