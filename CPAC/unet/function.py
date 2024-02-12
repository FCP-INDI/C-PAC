from click import BadParameter


class MyParser(BadParameter):
    def error(self, message):
        import sys
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        self.exit(2)


def write_nifti(data, aff, shape, out_path):
    import nibabel as nib

    data=data[0:shape[0],0:shape[1],0:shape[2]]
    img=nib.Nifti1Image(data, aff)
    img.to_filename(out_path)

def estimate_dice(gt_msk, prt_msk):
    intersection=gt_msk*prt_msk
    dice=2*float(intersection.sum())/float(gt_msk.sum()+prt_msk.sum())
    
    return dice

def extract_large_comp(prt_msk):
    
    import scipy.ndimage as snd
    import numpy as np

    labs, num_lab=snd.label(prt_msk)
    c_size=np.bincount(labs.reshape(-1))
    c_size[0]=0
    max_ind=c_size.argmax()
    prt_msk=labs==max_ind

    return prt_msk

def fill_holes(prt_msk):

    from CPAC.unet.function import extract_large_comp

    non_prt_msk=prt_msk==0
    non_prt_msk=extract_large_comp(non_prt_msk)
    prt_msk_filled=non_prt_msk==0

    return prt_msk_filled

def erosion_dilation(prt_msk, iterations=1):

    import scipy.ndimage as snd
    from CPAC.unet.function import extract_large_comp

    # Erosion
    structure=snd.generate_binary_structure(3, 1)
    prt_msk_eroded=snd.binary_erosion(prt_msk, structure=structure, iterations=iterations).astype(prt_msk.dtype)

    # Extract Largest Component
    prt_msk_eroded=extract_large_comp(prt_msk_eroded)

    # Dilation
    prt_msk_dilated=snd.binary_dilation(prt_msk_eroded, structure=structure, iterations=iterations).astype(prt_msk.dtype)

    return prt_msk_dilated

def predict_volumes(model_path, rimg_in=None, cimg_in=None, bmsk_in=None, suffix="unet_pre_mask", 
        ed_iter=0, save_dice=False, save_nii=True, nii_outdir=None, verbose=False, 
        rescale_dim=256, num_slice=3):

    import torch
    import torch.nn as nn
    import numpy as np
    from torch.autograd import Variable
    from CPAC.unet.function import extract_large_comp, estimate_dice, write_nifti, fill_holes, erosion_dilation
    from CPAC.unet.model import UNet2d
    from CPAC.unet.dataset import VolumeDataset, BlockDataset
    from torch.utils.data import DataLoader
    import os, sys
    import nibabel as nib
    import pickle

    train_model = UNet2d(dim_in=3, num_conv_block=5, kernel_root=16)
    checkpoint = torch.load(model_path, map_location={'cuda:0':'cpu'})
    train_model.load_state_dict(checkpoint['state_dict'])
    model = nn.Sequential(train_model, nn.Softmax2d())

    use_gpu=torch.cuda.is_available()
    model_on_gpu=next(model.parameters()).is_cuda
    use_bn=True
    if use_gpu:
        if not model_on_gpu:
            model.cuda()
    else:
        if model_on_gpu:
            model.cpu()

    NoneType=type(None)
    if isinstance(rimg_in, NoneType) and isinstance(cimg_in, NoneType):
        print("Input rimg_in or cimg_in")
        sys.exit(1)

    if save_dice:
        dice_dict=dict()
    
    volume_dataset=VolumeDataset(rimg_in=rimg_in, cimg_in=cimg_in, bmsk_in=bmsk_in)
    volume_loader=DataLoader(dataset=volume_dataset, batch_size=1)
    
    for idx, vol in enumerate(volume_loader):
        if len(vol)==1: # just img
            ptype=1 # Predict
            cimg=vol
            bmsk=None
            block_dataset=BlockDataset(rimg=cimg, bfld=None, bmsk=None, num_slice=num_slice, rescale_dim=rescale_dim)
        elif len(vol)==2: # img & msk
            ptype=2 # image test
            cimg=vol[0]
            bmsk=vol[1]
            block_dataset=BlockDataset(rimg=cimg, bfld=None, bmsk=bmsk, num_slice=num_slice, rescale_dim=rescale_dim)
        elif len(vol==3): # img bias_field & msk
            ptype=3 # image bias correction test
            cimg=vol[0]
            bfld=vol[1]
            bmsk=vol[2]
            block_dataset=BlockDataset(rimg=cimg, bfld=bfld, bmsk=bmsk, num_slice=num_slice, rescale_dim=rescale_dim)
        else:
            print("Invalid Volume Dataset!")
            sys.exit(2)
        
        rescale_shape=block_dataset.get_rescale_shape()
        raw_shape=block_dataset.get_raw_shape()
        
        for od in range(3):
            backard_ind=np.arange(3)
            backard_ind=np.insert(np.delete(backard_ind, 0), od, 0)

            block_data, slice_list, slice_weight=block_dataset.get_one_directory(axis=od)
            pr_bmsk=torch.zeros([len(slice_weight), rescale_dim, rescale_dim])
            if use_gpu:
                pr_bmsk=pr_bmsk.cuda()
            for (i, ind) in enumerate(slice_list):
                if ptype==1:
                    rimg_blk=block_data[i]
                    if use_gpu:
                        rimg_blk=rimg_blk.cuda()
                elif ptype==2:
                    rimg_blk, bmsk_blk=block_data[i]
                    if use_gpu:
                        rimg_blk=rimg_blk.cuda()
                        bmsk_blk=bmsk_blk.cuda()
                else:
                    rimg_blk, bfld_blk, bmsk_blk=block_data[i]
                    if use_gpu:
                        rimg_blk=rimg_blk.cuda()
                        bfld_blk=bfld_blk.cuda()
                        bmsk_blk=bmsk_blk.cuda()
                pr_bmsk_blk=model(torch.unsqueeze(Variable(rimg_blk), 0))
                pr_bmsk[ind[1], :, :]=pr_bmsk_blk.data[0][1, :, :]
            
            if use_gpu:
                pr_bmsk=pr_bmsk.cpu()
            
            pr_bmsk=pr_bmsk.permute(backard_ind[0], backard_ind[1], backard_ind[2])
            pr_bmsk=pr_bmsk[:rescale_shape[0], :rescale_shape[1], :rescale_shape[2]]
            uns_pr_bmsk=torch.unsqueeze(pr_bmsk, 0)
            uns_pr_bmsk=torch.unsqueeze(uns_pr_bmsk, 0)
            uns_pr_bmsk=nn.functional.interpolate(uns_pr_bmsk, size=raw_shape, mode="trilinear", align_corners=False)
            pr_bmsk=torch.squeeze(uns_pr_bmsk)

            if od==0:
                pr_3_bmsk=torch.unsqueeze(pr_bmsk, 3)
            else:
                pr_3_bmsk=torch.cat((pr_3_bmsk, torch.unsqueeze(pr_bmsk, 3)), dim=3)
        
        pr_bmsk=pr_3_bmsk.mean(dim=3)
        
        pr_bmsk=pr_bmsk.numpy()
        pr_bmsk_final=extract_large_comp(pr_bmsk>0.5)
        pr_bmsk_final=fill_holes(pr_bmsk_final)
        if ed_iter>0:
            pr_bmsk_final=erosion_dilation(pr_bmsk_final, iterations=ed_iter)
        
        if isinstance(bmsk, torch.Tensor):
            bmsk=bmsk.data[0].numpy()
            dice=estimate_dice(bmsk, pr_bmsk_final)
            if verbose:
                print(dice)

        t1w_nii=volume_dataset.getCurCimgNii()
        t1w_path=t1w_nii.get_filename()
        t1w_dir, t1w_file=os.path.split(t1w_path)
        t1w_name=os.path.splitext(t1w_file)[0]
        t1w_name=os.path.splitext(t1w_name)[0]

        if save_nii:
            t1w_aff=t1w_nii.affine
            t1w_shape=t1w_nii.shape

            if isinstance(nii_outdir, NoneType):
                nii_outdir = os.getcwd()

            out_path=os.path.join(nii_outdir, t1w_name+"_"+suffix+".nii.gz")
            write_nifti(np.array(pr_bmsk_final, dtype=np.float32), t1w_aff, t1w_shape, out_path)

        if save_dice:
            dice_dict[t1w_name]=dice

    if save_dice:
        return dice_dict
    
    # return output mask
    return out_path

