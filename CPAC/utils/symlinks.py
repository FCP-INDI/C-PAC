import os
import errno
from collections import defaultdict

output_renamings = {
    'anatomical_brain': 'anat',
    'anatomical_brain_mask': 'anat',
    'qc': 'qc',
    'anatomical_skull_leaf': 'anat',
    'anatomical_to_mni_linear_xfm': 'anat',
    'mni_to_anatomical_linear_xfm': 'anat',
    'mni_to_anatomical_nonlinear_xfm': 'anat',
    'anatomical_to_mni_nonlinear_xfm': 'anat',
    'anatomical_gm_mask': 'anat',
    'anatomical_csf_mask': 'anat',
    'anatomical_wm_mask': 'anat',
    'ants_initial_xfm': 'anat',
    'ants_rigid_xfm': 'anat',
    'ants_affine_xfm': 'anat',
    'mean_functional': 'func',
    'functional_preprocessed_mask': 'func',
    'functional_to_spatial_map': 'func',
    'functional_mask_to_spatial_map': 'func',
    'fmap_phase_diff': 'func',
    'fmap_magnitude': 'func',
    'functional_distortion_corrected': 'func',
    'despiked_fieldmap': 'func',
    'prepared_fieldmap_map': 'func',
    'fieldmap_mask': 'func',
    'slice_time_corrected': 'func',
    'slice_timing_corrected': 'func',
    'movement_parameters': 'parameters',
    'max_displacement': 'parameters',
    'xform_matrix': 'parameters',
    'output_means': 'parameters',
    'functional_preprocessed': 'func',
    'functional_brain_mask': 'func',
    'motion_correct': 'func',
    'motion_correct_smooth': 'func',
    'motion_correct_to_standard': 'func',
    'motion_correct_to_standard_smooth': 'func',
    'mean_functional_in_anat': 'func',
    'coordinate_transformation': 'func',
    'raw_functional': 'func',
    'selected_func_volume': 'func',
    'anatomical_wm_edge': 'registration',
    'anatomical_to_functional_xfm': 'registration',
    'inverse_anatomical_to_functional_xfm': 'registration',
    'functional_gm_mask': 'segmentation',
    'functional_wm_mask': 'segmentation',
    'functional_csf_mask': 'segmentation',
    'frame_wise_displacement_power': 'parameters',
    'frame_wise_displacement_jenkinson': 'parameters',
    'functional_nuisance_residuals': 'func',
    'functional_nuisance_regressors': 'func',
    'power_spectrum_distribution': 'alff',
    'functional_freq_filtered': 'func',
    'scrubbing_movement_parameters': 'parameters',
    'despiking_frames_included': 'parameters',
    'despiking_frames_excluded': 'parameters',
    'scrubbing_frames_included': 'parameters',
    'scrubbing_frames_excluded': 'parameters',
    'motion_params': 'parameters',
    'power_params': 'parameters',
    'scrubbed_preprocessed': 'func',
    'functional_to_standard': 'func',
    'functional_brain_mask_to_standard': 'func',
    'mean_functional_to_standard': 'func',
    'functional_to_anat_linear_xfm': 'registration',
    'functional_to_mni_linear_xfm': 'registration',
    'mni_to_functional_linear_xfm': 'registration',
    'ants_symmetric_initial_xfm': 'registration',
    'ants_symmetric_rigid_xfm': 'registration',
    'ants_symmetric_affine_xfm': 'registration',
    'anatomical_to_symmetric_mni_nonlinear_xfm': 'registration',
    'symmetric_mni_to_anatomical_nonlinear_xfm': 'registration',
    'symmetric_mni_to_anatomical_linear_xfm': 'registration',
    'anat_to_symmetric_mni_ants_composite_xfm': 'registration',
    'symmetric_anatomical_to_standard': 'registration',
    'anatomical_to_symmetric_mni_linear_xfm': 'registration',
    'anatomical_to_standard': 'anat',
    'leaf_node_to_standard': 'func',
    'vmhc_raw_score': 'vmhc',
    'vmhc_fisher_zstd': 'vmhc',
    'vmhc_fisher_zstd_zstat_map': 'vmhc',
    'alff': 'alff',
    'falff': 'alff',
    'alff_smooth': 'alff',
    'falff_smooth': 'alff',
    'alff_to_standard': 'alff',
    'falff_to_standard': 'alff',
    'alff_to_standard_smooth': 'alff',
    'falff_to_standard_smooth': 'alff',
    'alff_to_standard_zstd': 'alff',
    'falff_to_standard_zstd': 'alff',
    'alff_to_standard_smooth_zstd': 'alff',
    'falff_to_standard_smooth_zstd': 'alff',
    'alff_to_standard_zstd_smooth': 'alff',
    'falff_to_standard_zstd_smooth': 'alff',
    'reho': 'reho',
    'reho_smooth': 'reho',
    'reho_to_standard': 'reho',
    'reho_to_standard_smooth': 'reho',
    'reho_to_standard_zstd': 'reho',
    'reho_to_standard_smooth_zstd': 'reho',
    'reho_to_standard_zstd_smooth': 'reho',
    'voxel_timeseries': 'timeseries',
    'roi_timeseries': 'timeseries',
    'roi_timeseries_for_SCA': 'timeseries',
    'roi_timeseries_for_SCA_multreg': 'timeseries',
    'sca_roi_files': 'sca_roi',
    'sca_roi_files_smooth': 'sca_roi',
    'sca_roi_files_to_standard': 'sca_roi',
    'sca_roi_files_to_standard_smooth': 'sca_roi',
    'sca_roi_files_to_standard_fisher_zstd': 'sca_roi',
    'sca_roi_files_to_standard_smooth_fisher_zstd': 'sca_roi',
    'sca_roi_files_to_standard_fisher_zstd_smooth': 'sca_roi',
    'bbregister_registration': 'surface_registration',
    'left_hemisphere_surface': 'surface_registration',
    'right_hemisphere_surface': 'surface_registration',
    'vertices_timeseries': 'timeseries',
    'centrality': 'centrality',
    'centrality_smooth': 'centrality',
    'centrality_zstd': 'centrality',
    'centrality_smooth_zstd': 'centrality',
    'centrality_zstd_smooth': 'centrality',
    'centrality_graphs': 'centrality',
    'seg_probability_maps': 'anat',
    'seg_mixeltype': 'anat',
    'seg_partial_volume_map': 'anat',
    'seg_partial_volume_files': 'anat',
    'spatial_map_timeseries': 'timeseries',
    'spatial_map_timeseries_for_DR': 'timeseries',
    'dr_tempreg_maps_files': 'spatial_regression',
    'dr_tempreg_maps_files_smooth': 'spatial_regression',
    'dr_tempreg_maps_zstat_files': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_smooth': 'spatial_regression',
    'dr_tempreg_maps_files_to_standard': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_to_standard': 'spatial_regression',
    'dr_tempreg_maps_files_to_standard_smooth': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_to_standard_smooth': 'spatial_regression',
    'sca_tempreg_maps_files': 'sca_roi',
    'sca_tempreg_maps_files_smooth': 'sca_roi',
    'sca_tempreg_maps_zstat_files': 'sca_roi',
    'sca_tempreg_maps_zstat_files_smooth': 'sca_roi',
}

fork_ids = {
    'compcor': 'nuis',
    'selector': 'nuis',
}

def group_files_in_strategies(output_dir, paths):

    strategies = defaultdict(set)

    for path in paths:
        pieces = path.replace(output_dir, '').split(os.sep)
        related_strategy = tuple([
            p for p in pieces if any(
                p.startswith('_' + pattern) for pattern in fork_ids.keys()
            )
        ])
        strategies[related_strategy].add(path)

    # Get every file that is not affected by a strategy and add it to every other strategy
    strategies_keys = set(strategies.keys()) - set(tuple())
    for paths_without_strategy in strategies[tuple()]:
        for strategy in strategies_keys:
            strategies[strategy].add(paths_without_strategy)

    strategies = dict(strategies)

    # Add every file from primary strategy into derived strategies
    for strategy in strategies_keys:
        strategy_files = strategies[strategy]

        has_derived_strategy = False

        for specialized_strategy in strategies_keys:

            if specialized_strategy not in strategies:
                continue

            # specialized_strategy is derived from strategy
            if len(specialized_strategy) > len(strategy) and \
                strategy == specialized_strategy[0:len(strategy)]:

                strategies[specialized_strategy].update(strategy_files)
                has_derived_strategy = True

        if has_derived_strategy:
            del strategies[strategy]

    if tuple() in strategies:
        del strategies[tuple()]

    return strategies


def compile_strategy_name(strategy):
    name = []
    for s in strategy:
        fork = s[1:]
        for id in fork_ids:
            if fork.startswith(id):
                fork_name = fork_ids[id]
                name += [s.replace('_' + id, fork_name)]
                break
    return '__'.join(name)


def create_paths_to_symlinks(
    output_dir,
    pipeline_id,
    subject_id,
    paths
):

    if len(paths) == 0:
        return {}

    # Mapping of paths to symlinks
    symlinks = {}

    strategies = group_files_in_strategies(output_dir, paths)

    for strategy, strategy_paths in strategies.items():

        strategy_name = compile_strategy_name(strategy)

        paths_symlinks = []

        for path in strategy_paths:

            path_parts = path.replace(output_dir, '').strip(os.sep).split(os.sep)

            # If file is relative to a scan, use the tag
            # Otherwise, set to 'scan'
            scan_name = [p for p in path_parts if p.startswith('_scan_')]
            if not scan_name:
                scan_name = 'scan'
            else:
                scan_name = scan_name[0]
                path_parts = [p for p in path_parts if p != scan_name]
                scan_name = scan_name.strip('_')

            # Treat QC differently
            if path_parts[2] == 'qc':
                output_name = path_parts[3]
                output_group = 'qc'
                output_specs = path_parts[4:-1]
            else:
                output_name = path_parts[2]
                output_group = output_renamings[output_name]
                output_specs = path_parts[3:-1]
            # output_specs are parameters e.g. bandpass, fwhm, and mask

            # Remove strategies info from path, since it is included
            # on the outer path
            output_specs = [
                sp for sp in output_specs
                if not any(
                    sp.startswith('_' + pattern)
                    for pattern in fork_ids.keys()
                )
            ]

            # Remove leading or trailing underscores
            output_specs = [sp.strip('_') for sp in output_specs]

            # Get the file extension
            output_file_ext = ''
            output_file_parts = path_parts[-1].split('.')
            if output_file_parts[-2:] == ['nii', 'gz']:
                output_file_ext = '.nii.gz'
            else:
                output_file_ext = '.' + output_file_parts[-1]

            paths_symlinks.append({
                'path': path,
                'scan': scan_name,
                'group': output_group,
                'specs': tuple(output_specs),
                'name': output_name,
                'file': path_parts[-1],
                'ext': output_file_ext,
            })

        for i, path_symlink in enumerate(paths_symlinks):

            original_path = path_symlink['path']

            remaining_paths_symlinks = paths_symlinks[0:i] + paths_symlinks[i+1:]

            group = path_symlink['group']

            # Check if we need to keep the same basename or if we can give it
            # a better name, in order to avoid link collision
            # Outputs like from segmentation cant change its name since they
            # have several files (e.g. seg0, seg1, seg2)
            if any(
                rps['group'] == group and
                rps['scan'] == path_symlink['scan'] and
                rps['name'] == path_symlink['name'] and
                rps['specs'] == path_symlink['specs']
                for rps in remaining_paths_symlinks
            ):
                symlink_basename = path_symlink['file']

                # if pretty name cant be used, append it to group name
                group = os.path.join(group, path_symlink['name'])
            else:
                symlink_basename = path_symlink['name'] + path_symlink['ext']

            sym_path = os.path.join(*[
                pipeline_id,
                strategy_name,
                path_symlink['scan'],
                group,
            ] + list(path_symlink['specs']) + [symlink_basename])

            symlinks[original_path] = sym_path

    values = list(symlinks.values())
    duplicates = set([x for x in values if values.count(x) > 1])
    if duplicates:
        raise Exception("Found duplicates: " + str(duplicates))

    return symlinks


def create_symlinks(
    output_dir,
    pipeline_id,
    subject_id,
    paths,
    relative=True,
    symlink_dir='sym_links'
):

    original_cwd = os.getcwd()

    try:
        os.chdir(output_dir)

        mapping = create_paths_to_symlinks(
            output_dir,
            pipeline_id,
            subject_id,
            paths
        )

        for path, symlink in mapping.items():

            relpath = path
            if relpath.startswith(output_dir):
                relpath = relpath[len(output_dir):].lstrip('/')

            symlink = os.path.join(symlink_dir, symlink)

            try:
                os.makedirs(os.path.dirname(symlink))
            except:
                pass

            if relative:
                backtrack = os.path.join(*(
                    ['..'] * (len(symlink.split('/')) - 1)
                ))
                path = os.path.join(backtrack, relpath)

            try:
                os.symlink(path, symlink)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    os.remove(symlink)
                    os.symlink(path, symlink)
                else:
                    raise e

    finally:
        os.chdir(original_cwd)
