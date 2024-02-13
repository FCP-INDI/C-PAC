# Copyright (C) 2021-2024  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Get Package ID.

Script to get GHCR ID string for a given owner + image tag

Usage: python get_package_id.py $OWNER $IMAGE_TAG $VERSION_TAG
"""
import os
import sys

import requests


def get_packages(owner, tag, api_token=None):
    """Collect GHCR packages for a given owner & tag.

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
        api_token = os.environ.get("GITHUB_TOKEN", "")

    def fetch(url):
        """Make API call and return response, given a URL.

        Parameters
        ----------
        url : str

        Returns
        -------
        dict or list
        """
        response = requests.get(
            url, headers={"Authorization": f"token {api_token}"}
        ).json()
        if (
            isinstance(response, dict)
            and response.get("message", "") == "Bad credentials"
        ):
            raise PermissionError(
                "\n".join(
                    [
                        ": ".join(
                            [
                                response["message"],
                                api_token if api_token else "[no token provided]",
                            ]
                        ),
                        "Either set GITHUB_TOKEN to a personal access token with "
                        "read.packages permissions or explicitly pass one as a fourth "
                        "positional argument:\n"
                        "`python get_package_id.py $OWNER $IMAGE_TAG "
                        "$VERSION_TAG $GITHUB_TOKEN`",
                    ]
                )
            )
        return response

    _packages = fetch(
        f"https://api.github.com/orgs/{owner}/packages/container/" f"{tag}/versions"
    )
    packages = []
    for _package in _packages:
        if _package.get("message", "Not Found") == "Not Found":
            packages += fetch(
                f"https://api.github.com/users/{owner}/packages/container/"
                f"{tag}/versions"
            )
    return packages


def id_from_tag(owner, image, tag, api_token=None):
    """Return a package ID given an image version tag.

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
    versions = [
        image["id"]
        for image in packages
        if tag in image.get("metadata", {}).get("container", {}).get("tags", [])
    ]
    if len(versions):
        return versions[0]
    msg = f"Image not found: ghcr.io/{owner}/{image}:{tag}"
    raise LookupError(msg)


if __name__ == "__main__":
    if len(sys.argv) == 4:  # noqa: PLR2004
        print(id_from_tag(*sys.argv[1:]))  # noqa: T201
    else:
        print(__doc__)  # noqa: T201
