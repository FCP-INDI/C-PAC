import os
import random
import yaml

from sys import argv
from warnings import warn

platforms = {
    'docker': {
        'bind': '-v',
        'set_home': '--workdir=',
        'image': ':'.join([
            'fcpindi/c-pac',
            os.environ.get('CIRCLE_BRANCH', 'latest')
        ])
    },
    'singularity': {
        'bind': '-B',
        'set_home': '-H ',
        'image': 'C-PAC-CI.simg'
    }
}


def get_random_subject(species='human'):
    """
    Function to get a random data config file and subject for a given species.

    Note: only human data are configured at the moment

    Parameters
    ----------
    species: str
        'human', 'rodent', or 'nhp'

    Returns
    -------
    data_config_file: str

    participant_ndx: int
    """
    if species == 'human':
        data_config_file = 'test_configs/data-test_4-projects_5-subjects.yml'
    else:
        raise NotImplementedError(
            f'Data configurations not yet set for random test of {species}'
        )
    with open(data_config_file, 'r') as data_config:
        subject_list = yaml.safe_load(data_config)
    return (data_config_file, random.randrange(len(subject_list)))


def get_random_test_run_command(platform, image=None):
    """
    Function to choose a random preconfig, an appropriate subject, and
    return a shell command to run C-PAC on that subject/config in the
    given platform image.

    Parameters
    ----------
    platform: str

    image: str

    Returns
    -------
    command: str
    """
    if image is None:
        image = platforms[platform]['image']

    # collect preconfigs
    all_configs = {
        'default', 
        *{
            config[16:-4] for config in os.listdir(
                'CPAC/resources/configs'
            ) if config.startswith('pipeline_config')
        }
    }

    # choose a random preconfig
    random_config = random.choice(list(all_configs))
    config_string = '' if (
        random_config == 'default'
    ) else f'--preconfig {random_config}'

    # determine appropriate species
    if random_config in {'nhp-macaque', 'monkey'}:
        data_species = 'nhp'
    elif random_config == 'rodent':
        data_species = 'rodent'
    else:
        data_species = 'human'

    try:
        data_config_file, participant_ndx = get_random_subject(data_species)
        command = ' '.join([
            f'{platform} exec',
            f'{platforms[platform]["set_home"]}/home/circleci/project',
            platforms[platform]['bind'],
            'CPAC/resources/configs/test_configs:/test_configs',
            ' '.join([
                '-e COVERAGE_FILE=.coverage.docker-test-run',
                image
            ]) if platform == 'docker' else image,
            'coverage run /code/dev/docker_data/run.py',
            '/home/circleci/project',
            '/home/circleci/project/outputs participant',
            f'--save_working_dir --data_config_file {data_config_file}',
            f'{config_string} --n_cpus 1 --mem_gb 12'.lstrip(),
            f'--participant_ndx {participant_ndx}'
        ])
        if platform == 'singularity':
            command = ' '.join([
                'SINGULARITYENV_COVERAGE_FILE=.coverage.singularity-test-run',
                command
            ])
    except NotImplementedError as nie:
        # pass error along to user as warning, but don't fail
        warn(nie, Warning)
        command = (
            f'echo "{nie} which we need for \'--preconfig {random_config}\'"'
        )

    return command
    

if __name__ == '__main__':
    print(get_random_test_run_command(argv[1]))
