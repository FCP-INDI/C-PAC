import os
import random
import stat
import yaml

from warnings import warn


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
        data_config_file = (
            'CPAC/resources/configs/test_configs/'
            'data-test_4-projects_5-subjects.yml'
        )
    else:
        raise NotImplementedError(
            f'Data configurations not yet set for random test of {species}'
        )
    with open(data_config_file, 'r') as data_config:
        subject_list = yaml.safe_load(data_config)
    return (data_config_file, random.randrange(len(subject_list)))


def get_random_test_run_command():
    """
    Function to choose a random preconfig, an appropriate subject, and
    return a string command to pass to coverage_run.sh.

    Parameters
    ----------
    None

    Returns
    -------
    command: str
    """
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
            'python -m coverage run /code/dev/docker_data/run.py',
            '/home/circleci/project',
            '/home/circleci/project/outputs participant',
            f'--save_working_dir --data_config_file {data_config_file}',
            f'{config_string} --n_cpus 1 --mem_gb 12'.lstrip(),
            f'--participant_ndx {participant_ndx}'
        ])
    except NotImplementedError as nie:
        # pass error along to user as warning, but don't fail
        warn(nie, Warning)
        command = (
            f'echo "{nie}, which we need for \'--preconfig {random_config}\'"'
        )

    return command


if __name__ == '__main__':
    fp = os.path.join(os.path.dirname(__file__), 'run_command.sh')
    run_string = get_random_test_run_command()
    with open(fp, 'w') as run_command:
        run_command.write('\n'.join([
            '#!/bin/bash\n',
            '# install testing requirements',
            'pip install -r /code/dev/circleci_data/requirements.txt\n\n'
        ]))
        if run_string[:4] != 'echo':
            run_command.write('\n'.join([
                '# run one participant with coverage',
                run_string,
                ''
            ]))
        else:
            with open('run_warning.py', 'w') as warning_script:
                warning_script.write(f'print({run_string[5:]})')
            run_command.write('python -m coverage run run_warning.py')
    # ↓ chmod +x ↓
    os.chmod(fp, os.stat(fp).st_mode | 0o0111)