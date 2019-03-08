import yaml
import tempfile

from CPAC.utils.yaml_template import create_yaml_from_template


def test_yaml_template():

    config_file = tempfile.mkstemp(suffix='test_yaml_template')[1]

    with open("./default_pipeline.yaml", "r") as f:
        config = yaml.load(f)

    new_config = create_yaml_from_template(config, "./default_pipeline.yaml")

    with open(config_file, "wb") as f:
        f.write(new_config)

    # TODO test cases
    # Test lists, dicts and list of dicts
    # Look for who is in flow style and who is not
    # Check for default values
    # Check for comments