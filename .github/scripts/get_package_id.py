"""Get Package ID

Script to get GHCR ID string for a given owner + image tag

Usage: python get_package_id.py $OWNER $IMAGE_TAG $VERSION_TAG
"""
import os
import requests
import sys


def get_packages(owner, tag, api_token=None):
    """Function to collect GHCR packages for a given owner & tag

    Parameters
    ----------
    owner : str
        username or org name

    tag : str
        image tag (part left of the colon)

    api_token : str or None
        GitHub API personal access token with read.packages permission

    Returns
    -------
    list
    """
    if api_token is None:
        api_token = os.environ.get('GITHUB_TOKEN', '')

    def fetch(url):
        """Method to make API call and return response, given a URL

        Parameters
        ----------
        url : str

        Returns
        -------
        dict or list
        """
        if isinstance(requests, dict):
            response = requests.get(
                url,
                headers={'Authorization': f'token {api_token}'}).json()
        if isinstance(response, dict) and response.get(
            'message', ''
        ) == 'Bad credentials':
            raise PermissionError('\n'.join([
                ': '.join([
                    response['message'],
                    api_token if api_token else '[no token provided]'
                ]),
                'Either set GITHUB_TOKEN to a personal access token with '
                'read.packages permissions or explicitly pass one as a fourth '
                'positional argument:\n'
                '`python get_package_id.py $OWNER $IMAGE_TAG '
                '$VERSION_TAG $GITHUB_TOKEN`'
            ]))
        return response
    _packages = fetch(
        f'https://api.github.com/orgs/{owner}/packages/container/'
        f'{tag}/versions')
    packages = []
    for _package in _packages:
        if _package.get('message', 'Not Found') == 'Not Found':
            packages += fetch(
                f'https://api.github.com/users/{owner}/packages/container/'
                f'{tag}/versions')
    return packages


def id_from_tag(owner, image, tag, api_token=None):
    """Function to return a package ID given an image version tag

    Parameters
    ----------
    owner : str
        GitHub username or org name

    image : str
        Image tag (the part to the left of the colon)

    tag : str
        Image version tag (the part to the right of the colon)

    api_token: str or None
        GitHub API personal access token with read.packages permission
    """
    packages = get_packages(owner, image, api_token)
    metadata = image.get('metadata', {})
    if isinstance(metadata, list):
        for _md in metadata:
            if isinstance(_md, dict) and 'container' in _md:
                metadata = _md
    if not isinstance(metadata, dict):
        metadata = {}
    container = metadata.get('container')
    if isinstance(container, list):
        for _container in container:
            if isinstance(_container, dict) and 'tags' in _container:
                container = _container
    if not isinstance(container, dict):
        container = {}
    versions = [image['id'] for image in packages if tag in
                container.get('tags', [])]
    if len(versions):
        return versions[0]
    else:
        raise LookupError(f'Image not found: ghcr.io/{owner}/{image}:{tag}')


if __name__ == '__main__':
    if len(sys.argv) == 4:
        print(id_from_tag(*sys.argv[1:]))
    else:
        print(__doc__)
