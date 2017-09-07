# CPAC/utils/bids_metadata.py
#
# Author: Chris Gorgolewski
# Contributions: Daniel Clark

'''
'''


def get_metadata_for_nifti(bids_root, rel_path):
    '''
    '''

    # Get metadata based on bids nifti path

    # Import packages
    import json
    # TODO: port this to versions of fs after 0.5.4
    from fs.opener import fsopendir

    # Init variables
    fs = fsopendir(bids_root)

    rel_path = fs._decode_path(rel_path)
    sidecarJSON = rel_path.replace(".nii.gz", ".json")

    pathComponents = sidecarJSON.split('/')
    filenameComponents = pathComponents[-1].split("_")
    sessionLevelComponentList = [];
    subjectLevelComponentList = [];
    topLevelComponentList = [];
    ses = None;
    sub = None;

    for filenameComponent in filenameComponents:
        if filenameComponent[:3] != "run":
            sessionLevelComponentList.append(filenameComponent)
            if filenameComponent[:3] == "ses":
                ses = filenameComponent
            else:
                subjectLevelComponentList.append(filenameComponent)
                if filenameComponent[:3] == "sub":
                    sub = filenameComponent
                else:
                    topLevelComponentList.append(filenameComponent)

    topLevelJSON = "/" + "_".join(topLevelComponentList);
    potentialJSONs = [topLevelJSON]

    subjectLevelJSON = "/" + sub + "/" + "_".join(subjectLevelComponentList)
    potentialJSONs.append(subjectLevelJSON)

    if ses:
        sessionLevelJSON = "/" + sub + "/" + ses + "/" + "_".join(sessionLevelComponentList)
        potentialJSONs.append(sessionLevelJSON)

    potentialJSONs.append(sidecarJSON)

    merged_param_dict = {};
    for json_file_path in potentialJSONs:
        if fs.exists(json_file_path):
            param_dict = json.load(fs.open(json_file_path, "r"))
            merged_param_dict.update(param_dict)

    return merged_param_dict

