import yaml
import tempfile
import yamlordereddictloader

from CPAC.utils.yaml_template import create_yaml_from_template


def test_yaml_template():

    config_file = tempfile.mkstemp(suffix='test_yaml_template')[1]

    # Create a new YAML configuration file based on the default_pipeline.yml file.
    config = yaml.safe_load(open('./data/default_pipeline.yml', 'r'))

    new_config = create_yaml_from_template(config, './data/default_pipeline.yml')

    with open(config_file, 'wb') as f:
        f.write(new_config)

    # Verify that the output has preserved blank lines, comments
    with open(config_file, 'r') as f:
        lines = f.readlines()
        # Assert first lines starts with a comment
        assert lines[0].startswith('# ')

        # Assert that there are blank lines
        assert lines[6].isspace()
        assert lines[7].isspace()

        # Assert that regressors configuration written in block style
        assert lines[669].startswith('Regressors:')
        assert lines[670].startswith('- Bandpass:')

        # Assert that other configurations that are lists of values are written in flow style
        assert lines[764].startswith('roiTSOutputs: [true, true]')
