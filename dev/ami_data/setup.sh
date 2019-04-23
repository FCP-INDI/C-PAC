#!/bin/bash

set -e

SETUP_SCRIPT=""

function cleanup {
    shutdown -h now
}

trap cleanup EXIT

apt-get update && apt-get install -y awscli git jq

INSTANCE_ID="`wget -qO- http://instance-data/latest/meta-data/instance-id`"
REGION="`wget -qO- http://instance-data/latest/meta-data/placement/availability-zone | sed -e 's:\([0-9][0-9]*\)[a-z]*\$:\\1:'`"
AMI_NAME="`aws ec2 --region=${REGION} describe-tags --filters "Name=resource-id,Values=$INSTANCE_ID" "Name=key,Values=ami_name" --output text | cut -f5`"
AMI_DESCRIPTION="`aws ec2 --region=${REGION} describe-tags --filters "Name=resource-id,Values=$INSTANCE_ID" "Name=key,Values=ami_description" --output text | cut -f5`"

# Running tag
until aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags "Key=hey,Value=hou" --output json
do
    sleep 2
    echo "Waiting for roles"
done


git clone https://github.com/FCP-INDI/C-PAC.git /opt/C-PAC

if [ -z "${SETUP_SCRIPT}" ]
then
    bash /opt/C-PAC/dev/ami_data/setup_cpac.sh
else
    echo ${SETUP_SCRIPT} | base64 --decode > /tmp/setup_cpac.sh
    bash /tmp/setup_cpac.sh
    rm /tmp/setup_cpac.sh
fi

if [ $? != 0 ]; then
    until aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags "Key=ami,Value=error" --output json
    do
        sleep 2
        echo "Waiting for error tag"
    done
    exit 1
fi

aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags Key=ami,Value=setup --output json

IMAGE_ID=`aws ec2 --region=${REGION} create-image --instance-id ${INSTANCE_ID} --name "${AMI_NAME}" --description "${AMI_DESCRIPTION}" --no-reboot --output json | jq -r '.ImageId'`

echo "Waiting image building"
until [ `aws ec2 --region=${REGION} describe-images --filters "Name=image-id,Values=${IMAGE_ID}" --output json | jq -r '.Images | .[].State'` = 'available' ]; do
    sleep 2
done

aws ec2 --region=${REGION} create-tags --resources ${INSTANCE_ID} --tags Key=ami,Value=done --output json
echo "Image build!"