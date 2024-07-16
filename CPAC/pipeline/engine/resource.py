
from bids2table import BIDSTable, BIDSFile, join_bids_path, parse_bids_entities
from dataclasses import dataclass
import pandas as pd

class Resource(BIDSFile):
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
        self.metadata = row['meta__json'] if isinstance(row['meta__json'], dict) else {}
        self.row = row
        #extra stuff
        self.ds = {k: v for k, v in row.items() if k.startswith('ds')}
        self.entity = {k: v for k, v in row.items() if k.startswith('ent')}
        self.finfo = {k: v for k, v in row.items() if k.startswith('finfo')}
        
        self.filename = self.finfo['finfo__file_path'].split("/")[-1]
        self.file_path = self.finfo['finfo__file_path']
        self.rel_path = f"sub-{self.entity['ent__sub']}"
        if self.entity['ent__ses'] != "None":
            self.rel_path += f"/ses-{self.entity['ent__ses']}"
        self.rel_path += f"/{self.entity['ent__datatype']}"

        self.category = self.entity['ent__suffix']

        super().__init__(
            self.ds['ds__dataset'],
            self.ds['ds__dataset_path'],
            self.finfo['finfo__file_path'],
            self.entity,
            self.metadata,           
        )

        self.entity_to_bids_key = {
             
            # These are the keys that are used in generating the Resource names.
            # The order of the keys does matter and affects the Resource names.
            # Some of the keys can me commented or uncommented to include or exclude them from the Resource names.
            # If the keys are duplicated then the values are appended as a list.

            #'ent__sub': 'sub',
            #'ent__ses': 'ses',
            'ent__sample': 'sample',
            #'ent__task': 'task', 
            #'ent__acq': 'acq', 
            'ent__ce': 'ce', 
            'ent__trc': 'trc', 
            'ent__stain': 'stain', 
            'ent__rec': 'rec', 
            'ent__dir': 'dir',  
            'ent__mod': 'mod', 
            'ent__echo': 'echo', 
            'ent__flip': 'flip', 
            'ent__inv': 'inv', 
            'ent__mt': 'mt', 
            'ent__part': 'part', 
            'ent__proc': 'proc', 
            'ent__hemi': 'hemi', 
            'ent__space': 'space', 
            'ent__split': 'split', 
            'ent__recording': 'recording', 
            'ent__chunk': 'chunk', 
            'ent__res': 'res', 
            'ent__den': 'den', 
            'ent__label': 'label',
            'ent__extra_entities': 'extra_entities', 
            'ent__desc': 'desc',  
            'ent__suffix': 'suffix',
            #'ent__ext': 'ext'
        }
        
        self.name = self.generate_resource_name()
        self.strats = {
            str(self.cpac_provenance) : self.file_path
        }
        for key, value in self.metadata.items():
            setattr(self, key, value)
    
    def generate_resource_name(self):
        bids_key_value_pairs = []
        for col, bids_key in self.entity_to_bids_key.items():
            if pd.notnull(self.row[col]) and self.row[col] != '':
                value = self.row[col] if pd.notna(self.row[col]) else {}
                if bids_key in ['suffix']:
                    bids_key_value_pairs.append(value)
                elif bids_key in ['extra_entities']:
                    for key, value in value.items():
                        bids_key_value_pairs.append(f"{key}-{value}")
                else:
                    bids_key_value_pair = f"{bids_key}-{value}"
                    bids_key_value_pairs.append(bids_key_value_pair)
        bids_style_combined = '_'.join(bids_key_value_pairs)
        return bids_style_combined

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
            shutil.copy(self.finfo['finfo__file_path'], path_to_write)
            # Write the JSON file only if the ext is .nii.gz
            if self.filename.endswith('.nii.gz'):
                json_path = os.path.join(path_to_write, f"{self.filename.replace('.nii.gz', '.json')}")
                with open(json_path, 'w') as f:
                    f.write(json.dumps(self.metadata, indent=4))
            return f"successfully written to {path_to_write}"
        except Exception as e:
            WFLOGGER.error(f"Error writing to disk: {e}")
            print(f"Error writing to disk: {e}")
        