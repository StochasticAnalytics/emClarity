## BH_fftShift &rarr; EMC_maskIndex

### Signature:
```matlab
INDEX = BH_fftShift(windowRadius, SIZE, useGPU, varargin)
% vs
INDEX = EMC_maskIndex(TYPE, SIZE, METHOD, OPTION)
```

### Noticeable difference:
1. ```TYPE = 'half2full'```: ```EMC_maskIndex``` can compute the fftshift and ifftshift mask (as ```BH_fftshift```), but it can also compute a mask used by ```EMC_irfftn``` to switch from a half (non-redundant) to a full grid (redundant). With the ```fourierTransformer``` and the CUDA ```cufftExecR2C```/```cufftExecC2R``` it is not useful, but as the ```fourierTransformer``` is not tested, I'd rather stay with this for now.
    ```matlab
    % Compute the R2C:
    shape = [128,128];
    img1 = rand(size);
    img1_fft = EMC_rfftn(img1);  % half grid

    % ... do some stuff on the half grid ...

    % Compute the C2R:
    img2 = EMC_irfftn(img1_fft, shape);

    % So this ...
    ifftn(fftn(img1));
    % ... is equivalent to this:
    EMC_irfftn(EMC_fftn(img1), shape);
    
    % note: 'half2full' is only for not-centered grids at the moment.
    %       If it is useful, I[TF]'ll add this option for centered grids.
    ```
    ```EMC_irfftn``` is calling ```EMC_maskIndex``` to compute the 'half2full' index wrap mask. To speed things up, ```EMC_irfftn``` is saving the mask in the heap (persistent). As such, one should clear the function after use to release the allocated memory (```clear EMC_irfftn```). See ```help EMC_irfftn``` for more information.

2. ```fftshift``` and ```ifftshift```:
```EMC_maskIndex``` fully supports fftshift and ifftshift index masks, for 2d/3d full/half grid. It is also *considerably* faster than ```BH_fftshift``` as it doesn't rely on sub2ind and ndgrids.
    ```matlab
    % fftshift
    idx = BH_fftShift(0, [128,128], 0);
    % - equivalent -
    idx = EMC_maskIndex('fftshift', [128,128], 'cpu', {});
    
    % ifftshift
    idx = BH_fftShift(0, [-128,-128], 0);
    % - equivalent -
    idx = EMC_maskIndex('ifftshift', [128,128], 'cpu', {});
    
    % half grid
    idx = BH_fftShift(0, [-128,-128], 0, 'halfgrid');
    % - equivalent -
    idx = EMC_maskIndex('ifftshift', [128,128], 'cpu', {'half', true});
    
    % note: I [TF] don't think BH_fftShift is working properly for odd dimensions 3d half grids!
    ```
  
3. **windowRadius** is currently not supported, but will be added under ```OPTION.radius```.