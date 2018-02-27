

def read_group_list_text_file(group_list_text_file):
    """Read in the group-level analysis participant-session list text file."""

    with open(group_list_text_file, "r") as f:
        group_list = f.readlines()

    # each item here includes both participant and session, and this also will
    # become the main ID column in the written design matrix CSV
    group_list = [str(x).rstrip("\n") for x in group_list if x != ""]

    return group_list


def read_pheno_csv_into_df(pheno_csv):
    """Read the phenotypic file CSV into a Pandas DataFrame."""

    import pandas as pd

    with open(pheno_csv, "r") as f:
        pheno_df = pd.read_csv(f)

    return pheno_df


def get_pheno_evs(pheno_csv, ev_selections):
    """Get the columns of phenotypic data matched with their participant-
    session labels as a Pandas DataFrame."""
    pass


def create_design_matrix_df(group_list, pheno_df=None,
                             ev_selections=None, pheno_sub_label=None):
    """Create the design matrix intended for group-level analysis via the FSL
    FLAME tool.

    This does NOT create the final .mat file that FSL FLAME takes in. This is
    an intermediary design matrix CSV meant for the user to review.

    If there is a phenotype CSV provided, this function will align the
    participant-session ID labels in the CPAC individual-level analysis output
    directory with the values listed in the phenotype file.
    """

    import pandas as pd

    # map the participant-session IDs to just participant IDs
    group_list_map = {}
    for part_ses in group_list:
        group_list_map[part_ses.split("_")[0]] = part_ses

    # initialize a Pandas DataFrame for the new design matrix
    design_df = pd.DataFrame()

    # initialize the rows (no columns yet!)
    design_df["Participant-Session"] = group_list

    if pheno_df:
        # if a phenotype CSV file is provided with the data

        # TODO: next- this aligning of the pheno and the design DF
        # TODO: use the group_map to line up the pheno_sub_label column in the
        # TODO: pheno DF with the design DF rows- then grab your desired EV
        # TODO: column(s) and merge it with the design DF basically

        # align the pheno's participant ID column with the group sublist text
        # file
        if pheno_sub_label:
            pass

        # add any selected EVs to the design matrix from the phenotype file
        for ev in ev_selections:
            pass

    return design_df


def create_contrasts_template_df(design_df):
    """Create the template Pandas DataFrame for the contrasts matrix CSV.

    The headers in the contrasts matrix needs to match the headers of the
    design matrix."""
    pass


def preset_single_group_avg(group_list, pheno_df=None, covariate=None,
                            pheno_sub_label=None):
    """Set up the design matrix CSV for running a single group average
    (one-sample T-test)."""

    ev_selections = None
    if pheno_df and covariate and pheno_sub_label:
        # if we're adding an additional covariate
        ev_selections = [covariate]

    design_df = create_design_matrix_df(group_list, pheno_df,
                                        ev_selections=ev_selections,
                                        pheno_sub_label=pheno_sub_label)

    design_df["Group Mean"] = 1

    return design_df


def write_dataframe_to_csv(design_df):
    """Write out a matrix Pandas DataFrame into a CSV file."""
    pass


def run(group_list_text_file):

    group_list = read_group_list_text_file(group_list_text_file)

