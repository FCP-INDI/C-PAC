#!/bin/bash

set -e

STORAGE=30
REGION=$1
REGION="us-east-1"

VERSION=`cat ./version`

AMI_NAME="C-PAC ${VERSION}"
AMI_DESCRIPTION="Configurable Pipeline for the Analysis of Connectomes - Version ${VERSION}"

function cleanup {
    echo "Cleaning up..."

    if [ -n "$INSTANCE_ID" ]; then
        echo "Cleaning up instance"
        aws --region "${REGION}" ec2 terminate-instances --instance-ids "${INSTANCE_ID}" > /dev/null
        aws ec2 wait instance-terminated --instance-ids "${INSTANCE_ID}"
    fi

    if [ -n "$IMAGE_STORAGE" ]; then
        echo "Cleaning up local storage"
        rm -rf "${IMAGE_STORAGE}"
    fi

    if [ -n "$INSTANCE_RANDOM_NAME" ]; then
        echo "Cleaning up security"
        aws --region "${REGION}" ec2 delete-key-pair --key-name "${INSTANCE_RANDOM_NAME}"
        aws ec2 delete-security-group --group-name ${INSTANCE_RANDOM_NAME}

        echo "Cleaning up roles"
        aws --region "${REGION}" iam remove-role-from-instance-profile --instance-profile-name ${INSTANCE_RANDOM_NAME}.profile --role-name ${INSTANCE_RANDOM_NAME}
        aws --region "${REGION}" iam delete-instance-profile --instance-profile-name ${INSTANCE_RANDOM_NAME}.profile
        aws --region "${REGION}" iam detach-role-policy --role-name ${INSTANCE_RANDOM_NAME} --policy-arn arn:aws:iam::aws:policy/job-function/SystemAdministrator
        aws --region "${REGION}" iam delete-role --role-name ${INSTANCE_RANDOM_NAME}
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
aws --region "${REGION}" iam create-role --role-name ${INSTANCE_RANDOM_NAME} --assume-role-policy-document file://${IMAGE_TEMP}.role > /dev/null
aws --region "${REGION}" iam attach-role-policy --role-name ${INSTANCE_RANDOM_NAME} --policy-arn arn:aws:iam::aws:policy/job-function/SystemAdministrator > /dev/null
aws --region "${REGION}" iam create-instance-profile --instance-profile-name ${INSTANCE_RANDOM_NAME}.profile > /dev/null
aws --region "${REGION}" iam add-role-to-instance-profile --role-name ${INSTANCE_RANDOM_NAME} --instance-profile-name ${INSTANCE_RANDOM_NAME}.profile > /dev/null

echo "Creating security"
KEYPAIR=`aws --region "${REGION}" ec2 create-key-pair --key-name ${INSTANCE_RANDOM_NAME}`
echo ${KEYPAIR} | jq -r '.KeyMaterial' > ${IMAGE_TEMP}.pem
chmod 600 ${IMAGE_TEMP}.pem

INSTANCE_GROUP=`aws ec2 create-security-group --group-name ${INSTANCE_RANDOM_NAME} --description ${INSTANCE_RANDOM_NAME} | jq -r '.GroupId'`
aws ec2 authorize-security-group-ingress --group-name ${INSTANCE_RANDOM_NAME} --protocol tcp --port 22 --cidr 0.0.0.0/0

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
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=cpac-ami}, {Key=ami_name,Value=${AMI_NAME}}, {Key=ami_description,Value='${AMI_DESCRIPTION}'}]" \
    --instance-initiated-shutdown-behavior terminate \
    --output json
`
INSTANCE_ID=`echo ${INSTANCE_DATA} | jq -r '.Instances | .[].InstanceId'`

echo "Waiting for roles"
until aws --region "${REGION}" ec2 associate-iam-instance-profile --iam-instance-profile Name=${INSTANCE_RANDOM_NAME}.profile --instance-id ${INSTANCE_ID} > /dev/null 2>&1
do
    sleep 2
done

STARTED_AT=`date +%s`
ELAPSED=0

echo "Waiting for running tag"
while [ "$ELAPSED" -lt "120" ]; do
    TAG=`aws  --region "${REGION}" ec2 describe-tags --filters "Name=resource-id,Values=${INSTANCE_ID}" "Name=key,Values=hey" --output json | jq -r '.Tags | .[].Value'`
    if [ "${TAG}" = "hou" ]; then
        break
    fi
    NOW=`date +%s`
    ELAPSED=`expr ${NOW} - ${STARTED_AT}`
    sleep 2
done

if [ "$ELAPSED" -gt "120" ]; then
    echo "Instance timeout. Try again later."
    exit
else
    echo "Instance will create the image by itself! Waiting (it may take a while)."

    INSTANCE_IP=`aws --region "${REGION}" ec2 describe-instances --instance-ids ${INSTANCE_ID} \
                  --query "Reservations[*].Instances[*].PublicIpAddress" \
                  --output=text`


    function get_instance_logs {
        echo "ssh -o ConnectTimeout=3 -o ConnectionAttempts=1 -o StrictHostKeyChecking=no -i ${IMAGE_TEMP}.pem ubuntu@${INSTANCE_IP} tail -f /var/log/cloud-init-output.log"
        while true; do
            ssh -o ConnectTimeout=3 -o ConnectionAttempts=1 -o StrictHostKeyChecking=no -i ${IMAGE_TEMP}.pem ubuntu@${INSTANCE_IP} tail -f /var/log/cloud-init-output.log
        done
    }

    get_instance_logs &

    echo "Waiting for image setup"
    while true
    do
        ELAPSED=`expr ${NOW} - ${STARTED_AT}`
        STATUS=`aws --region "${REGION}" ec2 describe-instance-status --instance-ids ${INSTANCE_ID} | jq -r '.InstanceStatuses | .[].InstanceState.Name'`
        if [ "${STATUS}" != "running" ]; then
            echo "The image creation process stopped."
            exit
        fi
        if (( $ELAPSED > 14400 )); then  # 4 hours
            echo "It should not take that much (more then 4h). Please take a look."
            read 
        fi
        sleep 5
    done
fi