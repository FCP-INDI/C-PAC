
def test_check_s3():

    import os
    from CPAC.utils.datasource import check_for_s3

    data = check_for_s3(
        file_path='s3://fcp-indi/resources/cpac/resources/rois_2mm.nii.gz',
        creds_path=None,
        dl_dir='/tmp',
        img_type='anat'
    )

    assert os.path.isfile(data)


def test_check_s3_node():

    import os
    from CPAC.utils.datasource import create_check_for_s3_node

    node = create_check_for_s3_node(
        'image',
        file_path='s3://fcp-indi/resources/cpac/resources/mask-thr50-3mm.nii.gz',
        creds_path=None,
        dl_dir='/tmp',
        img_type='anat'
    )

    res = node.run()
    assert os.path.isfile(res.outputs.local_path)

    node = create_check_for_s3_node(
        'image',
        file_path='s3://fcp-indi/resources/cpac/resources/mask-thr50-3mm.nii.gz'
    )

    res = node.run()
    assert os.path.isfile(res.outputs.local_path)
