"""
Function to assume role in boto3, returns a boto3 S3 client with role setup.
"""

import boto3

def assume_role(role, user, duration=3600):

    """
    Assume the role required to do large data request
    :param role: role name to assume (last part of RoleArn)
    :param user: AWS username to assume role from
    :param duration: integer number of seconds to assume the role for
    :return: boto3 S3 client object
    """

    # Query the operator for an MFA token for use in the AWS authorisation

    mfa_token = input('Assume Role MFA token code: ')

    # create an STS client object that represents a live connection to the STS service
    sts_client = boto3.client('sts')

    # Call the assume_role method of the STSConnection object and pass the role ARN and a role session name.
    assumed_role_object = sts_client.assume_role(
        RoleArn=role,
        SerialNumber=user,
        RoleSessionName="AssumeRoleSession",
        DurationSeconds=duration,
        TokenCode=mfa_token
    )

    # Get the temporary credentials that can be used to make subsequent API calls
    credentials = assumed_role_object['Credentials']

    # Use the temporary credentials that AssumeRole returns to make a connection to Amazon S3
    client = boto3.client(
        's3',
        region_name='ap-southeast-2',
        aws_access_key_id=credentials['AccessKeyId'],
        aws_secret_access_key=credentials['SecretAccessKey'],
        aws_session_token=credentials['SessionToken'],
    )

    return client
