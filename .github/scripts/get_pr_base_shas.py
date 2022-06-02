import os
from github import Github

print(' '.join([
    pr.base.sha for pr in Github(os.environ.get(
        'GITHUB_TOKEN'
    )).get_repo(os.environ.get(
        'GITHUB_REPOSITORY'
    )).get_commit(os.environ.get(
        'GITHUB_SHA'
    )).get_pulls()
]))
