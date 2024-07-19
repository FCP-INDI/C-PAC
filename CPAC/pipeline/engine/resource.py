
from bids2table import BIDSTable, BIDSFile, join_bids_path, parse_bids_entities
from dataclasses import dataclass
import pandas as pd

class Resource():
    row: dict
    CpacProvenance: tuple
    ds: dict
    entity: dict
    finfo: dict
    metadata: dict
    entity_to_bids_key: dict
    name: str
    rel_path: str

    def __init__(self, row, CpacProvenance):

        self.cpac_provenance = CpacProvenance
        self.metadata = {} # replace with >> row['json'] if isinstance(row['json'], dict) else {}
        self.row = row
        for key, value in self.row.items():
            setattr(self, key, value)
        
        self.filename = self.file_path.split("/")[-1]
        self.rel_path = f"sub-{self.sub}"
        if self.ses != "None":
            self.rel_path += f"/ses-{self.ses}"
        self.rel_path += f"/{self.datatype}"

        self.suffix = self.suffix
        
        self.name = self.filename.split(".")[0]
        self.strats = {
            str(self.cpac_provenance) : self.file_path
        }
        for key, value in self.metadata.items():
            setattr(self, key, value)

    def __repr__(self):
        exclude_list = ['CpacConfig', 'CpacConfigHash', 'CpacProvenance', 'metadata', 'cpac_provenance', 'ds', 'entity', 'finfo', 'row', 'filename', 'file_path', 'rel_path', 'entities', 'path', 'entity_to_bids_key', ]  # Add attribute names to exclude
        attributes = {attr: value for attr, value in self.__dict__.items() if attr not in exclude_list and value is not None}
        return f"{self.__class__.__name__}({attributes})"

    # write to disk
    def write_to_disk(self, path):
        import shutil
        try:
            path_to_write = os.path.join(path, self.rel_path)
            os.makedirs(path_to_write, exist_ok=True)
            # Copy the NIFTI file
            shutil.copy(self.finfo['file_path'], path_to_write)
            # Write the JSON file only if the ext is .nii.gz
            if self.filename.endswith('.nii.gz'):
                json_path = os.path.join(path_to_write, f"{self.filename.replace('.nii.gz', '.json')}")
                with open(json_path, 'w') as f:
                    f.write(json.dumps(self.metadata, indent=4))
            return f"successfully written to {path_to_write}"
        except Exception as e:
            WFLOGGER.error(f"Error writing to disk: {e}")
            print(f"Error writing to disk: {e}")
        