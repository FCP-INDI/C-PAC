import torch
import torch.nn as nn

from torch.autograd import Variable

def weigths_init(m):
    if isinstance(m, nn.Conv2d):
        nn.init.normal_(m.weight.data)
        nn.init.fill_(m.bias.data)

def Conv3dBlock(dim_in, dim_out,
        kernel_size=3, stride=1, padding=1,
        bias=True, use_bn=False):
    if use_bn:
        return nn.Sequential(
            nn.Conv3d(dim_in,  dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.BatchNorm3d(dim_out),
            nn.LeakyReLU(0.1),
            nn.Conv3d(dim_out, dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.BatchNorm3d(dim_out),
            nn.LeakyReLU(0.1)
            )
    else:
        return nn.Sequential(
            nn.Conv3d(dim_in,  dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.LeakyReLU(0.1),
            nn.Conv3d(dim_out, dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.LeakyReLU(0.1)
            )
    
def UpConv3dBlock(dim_in, dim_out, 
        kernel_size=4, stride=2, padding=1,
        bias=False):
    return nn.Sequential(
        nn.ConvTranspose3d(dim_in, dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
        nn.LeakyReLU(0.1)
        )

def Conv2dBlock(dim_in, dim_out,
        kernel_size=3, stride=1, padding=1,
        bias=True, use_bn=True):
    if use_bn:
        return nn.Sequential(
            nn.Conv2d(dim_in,  dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.BatchNorm2d(dim_out),
            nn.LeakyReLU(0.1),
            nn.Conv2d(dim_out, dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.BatchNorm2d(dim_out),
            nn.LeakyReLU(0.1)
            )
    else:
        return nn.Sequential(
            nn.Conv2d(dim_in,  dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.LeakyReLU(0.1),
            nn.Conv2d(dim_out, dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
            nn.LeakyReLU(0.1)
            )
    
def UpConv2dBlock(dim_in, dim_out, 
        kernel_size=4, stride=2, padding=1,
        bias=True):
    return nn.Sequential(
        nn.ConvTranspose2d(dim_in, dim_out, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias),
        nn.LeakyReLU(0.1)
        )

class UNet3d(nn.Module):
    def __init__(self, 
            dim_in=1, num_conv_block=2, kernel_root=8, 
            use_bn=False):
        super(UNet3d, self).__init__()
        self.layers=dict()
        self.num_conv_block=num_conv_block
        # Conv Layers
        for n in range(num_conv_block):
            if n==0:
                setattr(self, "conv%d" % (n+1), Conv3dBlock(dim_in, kernel_root, use_bn=use_bn))
            else:
                setattr(self, "conv%d" % (n+1), Conv3dBlock(kernel_root*(2**(n-1)), kernel_root*(2**n), use_bn=use_bn))

        # UpConv Layers
        for n in range(num_conv_block-1):
            i=num_conv_block-1-n
            setattr(self, "upconv%dto%d" % (i+1, i), UpConv3dBlock(kernel_root*(2**i), kernel_root*(2**(i-1))))
            setattr(self, "conv%dm" % (i), Conv3dBlock(kernel_root*(2**i), kernel_root*(2**(i-1))))    
        setattr(self, "max_pool", nn.MaxPool3d(2))
        setattr(self, "out_layer", nn.Conv3d(kernel_root, 2, 3, 1, 1))
        
        # Weight Initialization
        for m in self.modules(): 
            if isinstance(m, nn.Conv3d) or isinstance(m, nn.ConvTranspose3d):
                m.weight.data.normal_(0, 0.02) 
                if m.bias is not None: 
                    m.bias.data.zero_() 
            elif isinstance(m, nn.BatchNorm3d): 
                m.weight.data.normal_(1.0, 0.02)

    def forward(self, x):
        num_conv_block=self.num_conv_block
        conv_out=dict()
        for n in range(num_conv_block):
            if n==0:
                conv_out["conv%d" % (n+1)]=getattr(self, "conv%d" % (n+1))(x)
            else:
                conv_out["conv%d" % (n+1)]=getattr(self, "conv%d" % (n+1))(self.max_pool(conv_out["conv%d" % n])) 

        for n in range(num_conv_block-1):
            i=num_conv_block-1-n
            tmp=torch.cat(
                    (
                    getattr(self, "upconv%dto%d" % (i+1, i))(conv_out["conv%d" % (i+1)]),
                    conv_out["conv%d" % (i)]
                    ),
                    1
                )
            out=getattr(self, "conv%dm" % (i))(tmp)

        out=self.out_layer(out)
        if not self.training:
            softmax_layer=nn.Softmax(dim=1)
            out=softmax_layer(out)
        return out

class UNet2d(nn.Module):
    def __init__(self, 
            dim_in=6, num_conv_block=3, kernel_root=4, 
            use_bn=True):
        super(UNet2d, self).__init__()
        self.layers=dict()
        self.num_conv_block=num_conv_block
        # Conv Layers
        for n in range(num_conv_block):
            if n==0:
                setattr(self, "conv%d" % (n+1), Conv2dBlock(dim_in, kernel_root, use_bn=use_bn))
            else:
                setattr(self, "conv%d" % (n+1), Conv2dBlock(kernel_root*(2**(n-1)), kernel_root*(2**n), use_bn=use_bn))

        # UpConv Layers
        for n in range(num_conv_block-1):
            i=num_conv_block-1-n
            setattr(self, "upconv%dto%d" % (i+1, i), UpConv2dBlock(kernel_root*(2**i), kernel_root*(2**(i-1))))
            setattr(self, "conv%dm" % (i), Conv2dBlock(kernel_root*(2**i), kernel_root*(2**(i-1))))    
        setattr(self, "max_pool", nn.MaxPool2d(2))
        #setattr(self, "out_layer", nn.Sequential(nn.Conv2d(kernel_root, 2, 3, 1, 1), nn.Softmax2d()))
        setattr(self, "out_layer", nn.Conv2d(kernel_root, 2, 3, 1, 1))
        
        # Weight Initialization
        self.apply(self.weights_init)

    def weights_init(self, m):
        if isinstance(m, nn.Conv2d) or isinstance(m, nn.ConvTranspose2d):
            m.weight.data.normal_(0, 0.02) 
            if m.bias is not None: 
                m.bias.data.zero_() 
        elif isinstance(m, nn.BatchNorm2d): 
            m.weight.data.normal_(1.0, 0.02)

    def forward(self, x):
        num_conv_block=self.num_conv_block
        conv_out=dict()
        for n in range(num_conv_block):
            if n==0:
                conv_out["conv%d" % (n+1)]=getattr(self, "conv%d" % (n+1))(x)
            else:
                conv_out["conv%d" % (n+1)]=getattr(self, "conv%d" % (n+1))(self.max_pool(conv_out["conv%d" % n])) 

        for n in range(num_conv_block-1):
            i=num_conv_block-1-n
            if n==0:
                tmp=torch.cat(
                        (
                        getattr(self, "upconv%dto%d" % (i+1, i))(conv_out["conv%d" % (i+1)]),
                        conv_out["conv%d" % (i)]
                        ),
                        1
                    )
            else:
                tmp=torch.cat(
                        (
                        getattr(self, "upconv%dto%d" % (i+1, i))(out),#(conv_out["conv%d" % (i+1)]),
                        conv_out["conv%d" % (i)]
                        ),
                        1
                    )

            out=getattr(self, "conv%dm" % (i))(tmp)

        out=self.out_layer(out)
        return out

class MultiSliceBcUNet(nn.Module):
    def __init__(self, 
            num_slice=6, in_shape=256,
            num_conv_block=4, kernel_root=16, 
            use_bn=True):
        super(MultiSliceBcUNet, self).__init__()
        
        for i in range(num_slice):
            setattr(self, "slice%d" % (i+1),
                nn.Sequential(
                    UNet2d(dim_in=num_slice, num_conv_block=num_conv_block, kernel_root=kernel_root, use_bn=use_bn),
                    nn.Conv2d(2, 1, kernel_size=1, stride=1, padding=0),
                    nn.ReLU()
                    )
                )

        self.num_slice=num_slice

    def forward(self, x):
        for i in range(self.num_slice):
            pho=getattr(self, "slice%d" % (i+1))(x)
            if i==0:
                out=pho
            else:
                out=torch.cat(
                    ( 
                        out, 
                        pho
                    ),
                    1
                    )

        return out
    
    def freeze(self):
        for param in model.parameters():
            param.requires_grad=False
    
    def unfreeze(self):
        for param in model.parameters():
            param.requires_grad=True

class MultiSliceSsUNet(nn.Module):
    def __init__(self, 
            num_slice=6, in_shape=256,
            num_conv_block=5, kernel_root=16, 
            use_bn=True):
        super(MultiSliceSsUNet, self).__init__()
        
        for i in range(num_slice):
            setattr(self, "slice%d" % (i+1), 
                UNet2d(dim_in=num_slice, num_conv_block=num_conv_block, kernel_root=kernel_root, use_bn=use_bn)
                )

        self.num_slice=num_slice

    def forward(self, x):
        for i in range(self.num_slice):
            pho=torch.unsqueeze(getattr(self, "slice%d" % (i+1))(x), 2)
            if i==0:
                out=pho
            else:
                out=torch.cat(
                    ( 
                        out, 
                        pho
                    ),
                    2
                    )

        return out
    
    def freeze(self):
        for param in model.parameters():
            param.requires_grad=False
    
    def unfreeze(self):
        for param in model.parameters():
            param.requires_grad=True

class MultiSliceModel(nn.Module):
    def __init__(self, 
            num_slice=6, in_shape=256,
            bc_num_conv_block=3, bc_kernel_root=8, 
            ss_num_conv_block=4, ss_kernel_root=8, 
            use_bn=True):
        super(MultiSliceModel, self).__init__()
        
        self.BcUNet=MultiSliceBcUNet(num_slice=num_slice, in_shape=in_shape,
            num_conv_block=bc_num_conv_block, kernel_root=bc_kernel_root,
            use_bn=use_bn)
        self.SsUNet=MultiSliceSsUNet(num_slice=num_slice, in_shape=in_shape,
            num_conv_block=ss_num_conv_block, kernel_root=ss_kernel_root,
            use_bn=use_bn)

    def forward(self, x, model='forward_full'):
        if model=="forward_bc_part":
            b_field=self.BcUNet(x)
            out=b_field
        elif model=="forward_ss_part":
            b_msk=self.SsUNet(x)
            out=b_msk
        elif model=="forward_full":
            b_field=self.BcUNet(x)
            x=x*b_field
        
            out=self.SsUNet(x)

        return out
