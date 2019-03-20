#!/bin/bash

function cleanup {
    shutdown -h now
}

trap cleanup EXIT

apt-get update && apt-get install awscli git

INSTANCE_ID="`wget -qO- http://instance-data/latest/meta-data/instance-id`"
REGION="`wget -qO- http://instance-data/latest/meta-data/placement/availability-zone | sed -e 's:\([0-9][0-9]*\)[a-z]*\$:\\1:'`"
AMI_NAME="`aws ec2 --region=${REGION} describe-tags --filters "Name=resource-id,Values=$INSTANCE_ID" "Name=key,Values=ami_name" --output text | cut -f5`"
AMI_DESCRIPTION="`aws ec2 --region=${REGION} describe-tags --filters "Name=resource-id,Values=$INSTANCE_ID" "Name=key,Values=ami_description" --output text | cut -f5`"

until aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags Key=hey,Value=hou --output json
do
    sleep 2
    echo "Waiting for roles"
done


git clone https://github.com/FCP-INDI/C-PAC.git /opt/C-PAC

bash /opt/C-PAC/dev/ami_data/setup_cpac.sh

if [ $? != 0 ]; then
    aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags Key=ami,Value=error --output json
    exit $ERROR_CODE
fi

aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags Key=ami,Value=setup --output json

IMAGE_ID=`aws ec2 create-image --instance-id ${INSTANCE_ID} --name ${AMI_NAME} --description ${AMI_DESCRIPTION} --no-reboot --output json | jq -r '.ImageId'`

until [ `aws ec2 --region=${REGION} describe-images --filters "Name=image-id,Values=ami-04e8edf5a023838dd" --output json | jq -r '.Images | .[].State'` = 'available' ]; do
    sleep 2
done

aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags Key=ami,Value=done --output json