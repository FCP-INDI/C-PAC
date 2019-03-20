#!/bin/bash

STORAGE=30
REGION=$1
REGION="us-east-1"

VERSION=`cat ./version`
VERSION='TESTING'

AMI_NAME="C-PAC ${VERSION}"
AMI_DESCRIPTION="Configurable Pipeline for the Analysis of Connectomes - Version ${VERSION}"

function cleanup {
    if [ -n "$INSTANCE_KEYPAIR_NAME" ]; then
        aws --region "${REGION}" ec2 create-key-pair --key-name "${INSTANCE_KEYPAIR_NAME}"
    fi

    if [ -n "$INSTANCE" ]; then
        aws --region "${REGION}" ec2 terminate-instances "${INSTANCE_ID}"
    fi

    if [ -n "$IMAGE_STORAGE" ]; then
        rm -rf "${IMAGE_STORAGE}"
    fi
}

trap cleanup EXIT

IMAGE_DATA=`
    aws --region "${REGION}" ec2 describe-images \
    --owners 099720109477 \
    --filters \
        'Name=name,Values=ubuntu/images/hvm-ssd/ubuntu-*-18.04-amd64-server-*' \
        'Name=state,Values=available' \
    --output json | jq -r '.Images | sort_by(.CreationDate) | last(.[])'
`

IMAGE_ID=`echo ${IMAGE_DATA} | jq -r '.ImageId'`
IMAGE_TEMP=`mktemp`
INSTANCE_RANDOM_NAME=`basename ${IMAGE_TEMP}`

cat <<EOT > ${IMAGE_TEMP}.role
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "ec2.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOT

echo "Creating roles"
aws iam create-role --role-name ${INSTANCE_RANDOM_NAME} --assume-role-policy-document file://${IMAGE_TEMP}.role
aws iam attach-role-policy --role-name ${INSTANCE_RANDOM_NAME} --policy-arn arn:aws:iam::aws:policy/job-function/SystemAdministrator
aws iam create-instance-profile --instance-profile-name ${INSTANCE_RANDOM_NAME}.profile
aws iam add-role-to-instance-profile --role-name ${INSTANCE_RANDOM_NAME} --instance-profile-name ${INSTANCE_RANDOM_NAME}.profile

aws --region "${REGION}" ec2 create-key-pair --key-name ${INSTANCE_RANDOM_NAME}

INSTANCE_GROUP=`
    aws --region "${REGION}" ec2 describe-security-groups --group-names default --output json | \
    jq -r '.SecurityGroups | last(.[]).GroupId'
`

echo ${IMAGE_DATA} | jq -r ".BlockDeviceMappings | map(select(.Ebs)) | .[].Ebs.VolumeSize = ${STORAGE}" > ${IMAGE_TEMP}.storage

echo "Running instance"
INSTANCE_DATA=`
aws --region "${REGION}" ec2 run-instances \
    --image-id ${IMAGE_ID} \
    --count 1 \
    --security-group-ids ${INSTANCE_GROUP} \
    --instance-type t2.medium  \
    --key-name ${INSTANCE_RANDOM_NAME} \
    --user-data file://dev/ami_data/setup.sh \
    --block-device-mappings file://${IMAGE_TEMP}.storage \
    --tag-specifications "ResourceType=instance,Tags=[{Key=ami_name,Value=${AMI_NAME}}, {Key=ami_description,Value='${AMI_DESCRIPTION}'}]" \
    --instance-initiated-shutdown-behavior terminate \
    --output json
`
INSTANCE_ID=`echo ${INSTANCE_DATA} | jq -r '.Instances | .[].InstanceId'`

until aws ec2 associate-iam-instance-profile --iam-instance-profile Name=${INSTANCE_RANDOM_NAME}.profile --instance-id ${INSTANCE_ID} > /dev/null 2>&1
do
    sleep 2
    echo "Waiting for roles"
done

STARTED_AT=`date +%s`
ELAPSED=0
while (( ELAPSED < 120 )); do
    TAG=`aws ec2 describe-tags --filters "Name=resource-id,Values=${INSTANCE_ID}" "Name=key,Values=hey" --output json | jq -r '.Tags | .[].Value'`
    if [ "${TAG}" = "hou" ]; then
        break
    fi
    echo "Waiting for running tag"
    NOW=`date +%s`
    ELAPSED=`expr ${NOW} - ${STARTED_AT}`
    sleep 2
done

if (( ELAPSED >= 120 )); then
    echo "Instance timeout"
else
    echo "Instance will create the image by itself! Bye."
fi