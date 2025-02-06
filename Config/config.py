import yaml

def load_config(config_file="Config/config.yaml"):
    with open(config_file, "r", encoding="utf-8") as file:
        return yaml.safe_load(file)