## BH_mask3d &rarr; EMC_maskShape and EMC_maskReference
BH_mask3d was breakdown into 2 functions:
- coordinates/EMC_maskShape: compute a real-space mask of a given shape.
- coordinates/EMC_maskReference:  compute a real-space reference mask of a 3d/2d object.

### Signature:
```matlab
[MASK, volCOM] = BH_mask3d(SHAPE, SIZE, RADIUS, CENTER, varargin);
% vs
MASK = EMC_maskShape(SHAPE, SIZE, RADIUS, METHOD, OPTION);
[MASK, COM, FRACTION] = EMC_referenceMask(IMAGE, PIXEL, METHOD, OPTION);
```

#### EMC_maskShape
1. Input format - 2d and 3d case
    ```matlab
    %% 3d
    % this ...
    MASK = BH_mask3d('sphere', [100,100,100], [30,30,30], [0,0,0]);
    % is equivalent to this:
    MASK = EMC_maskShape('sphere', [100,100,100], [30,30,30], 'gpu', {});
    
    %% 2d
    % this ...
    MASK = BH_mask3d('sphere', [100,100], [30,30], [0,0], '2d');
    % is equivalent to this:
    MASK = EMC_maskShape('sphere', [100,100], [30,30], 'gpu', {});
    ```

2. **shifts** -vs- **CENTER**: These are equivalent parameters, but as every EMC functions, ```EMC_maskShape``` keeps the origin and shifts as optional parameters.
    ```matlab
    % this ...
    MASK = BH_mask3d('sphere', [100,100], [30,30], [-20,10], '2d');
    % is equivalent to this:
    MASK = EMC_maskShape('sphere', [100,100], [30,30], 'gpu', {'shift', [-20,10]});
    ```

3. Asymmetric restriction - 'sym': Restrict the mask to the first asymmetric unit.
    ```matlab
    % this ...
    MASK = BH_mask3d('sphere', [100,100], [30,30], [0,0], '2d', 3);  % C3 sym
    % is equivalent to this:
    MASK = EMC_maskShape('sphere', [100,100], [30,30], 'gpu', {'sym', 3});  % C3 sym
    
    % note: EMC_maskShape supports 'sym' for every shape, not only for cylinders.
    ```

4. **Pixel roll off** - **'kernel'**: One big difference between EMC_maskShape and BH_mask3d is the size of the roll off (at the edge of the shape, between 0 and 1). BH_mask3d uses a fixed 6 pixels roll, which is not accessible from outside the function, whereas EMC_maskShape uses a roll off of variable length (by default, 5% of min(SIZE), with minimum of 8 pixels). Moreover, this parameters is fully accessible with the 'kernel' optional parameters.
    ```matlab
     MASK = EMC_maskShape('sphere', [100,100], [30,30], 'gpu', {'kernel', false});  % no roll off
     MASK = EMC_maskShape('sphere', [100,100], [30,30], 'gpu', {'kernel', 0.1});  % 10% of min(SIZE)
     
     % convolve the 'binary' mask with your own kernel (must be row vector; cf. documentation)
     ownKernel = makeMyOwnSeparableKernel();  % 1xn (n is odd)
     MASK = EMC_maskShape('sphere', [100,100], [30,30], 'gpu', {'kernel', ownKernel});
    ```

#### EMC_maskReference
1. **Global variables**: BH_mask3d uses global variable to set some parameters (see below). These are now optional parameters, accessible with the OPTION cell|struct.
    ```
    bh_global_binary_mask_low_pass  -> 'lowpass'
    bh_global_binary_mask_threshold -> 'threshold'
    bh_global_vol_est_scaling       -> 'hydration_scaling'
    ```
