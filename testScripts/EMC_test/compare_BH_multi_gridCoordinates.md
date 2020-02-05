## BH_multi_gridCoordinates &rarr; EMC_coordVectors, EMC_coordGrids and EMC_coordTransform
BH_multi_gridCoordinates was breakdown into 3 functions:
- coordinates/EMC_coordVectors: compute the coordinate vectors.
- coordinates/EMC_coordGrids: compute coordinate grids of different coordinate systems (cartesian, cylindrical, etc.)
- coordinates/EMC_coordTransform: compute the (transformed) cartesian grids used for interpolation.

### Signature:
```matlab
[Gc1,Gc2,Gc3,g1,g2,g3] = BH_multi_gridCoordinates(SIZE, SYSTEM, METHOD, TRANSFORMATION, ...
                                                  flgFreqSpace, flgShiftOrigin, flgRad, ...
                                                  varargin);
% vs
[g1,g2,g3] = EMC_coordVectors(SIZE, METHOD, OPTION);
[Gc1,Gc2,Gc3,g1,g2,g3] = EMC_coordGrids(SYSTEM, SIZE, METHOD, OPTION);
[Gc1,Gc2,Gc3,g1,g2,g3] = EMC_coordTransform(SIZE, METHOD, OPTION, varargin)
% note: EMC_coordTransform is currently the only EMC function that has varargin in place.
```

### Not supported:
1. **flgShiftOrigin** vs **origin**: Some origin conventions are not supported and new ones are added:
    ```matlab
    % size = [6, 7];
    
    %% SIMILAR
    % flgShiftOrigin = 1 is equivalent to origin = 1
    vX = [-3 -2 -1  0  1  2];
    vY = [-3 -2 -1  0  1  2  3];
    
    % flgShiftOrigin = 0 is ALMOST equivalent to origin = -1
    % BH:
    vX = [0  1  2  3 -2 -1];  % 3 ; probably an error?
    vY = [0  1  2  3 -3 -2 -1];
    % EMC:
    vX = [0  1  2 -3 -2 -1];  % -3
    vY = [0  1  2  3 -3 -2 -1];
    
    %% NOT SUPPORTED
    % flgShiftOrigin = -1 is not supported (I [TF] don't know what they are used for)
    vX = [1 2 3 4 5 6];
    vY = [1 2 3 4 5 6 7];
    
    % flgShiftOrigin = -2 is not supported (I [TF] don't know what they are used for)
    vX = [4 5 6 1 2 3];
    vY = [5 6 7 1 2 3 4];

    % ADDITION CONVENTIONS
    % origin = 0
    vX = [-2.5 -1.5 0 1.5 2.5];
    vY = [-3 -2 -1  0  1  2  3];
    
    % origin = 2
    vX = [-2 -1  0  1  2  3];
    vY = [-3 -2 -1  0  1  2  3];
    ```
2. **Scaling**|**Mag**: EMC_coordVectors and EMC_coordGrids do NOT support magnification changes, but EMC_coordTransform does.

---
### Noticeable differences

#### EMC_coordVectors
1. **Notation** - **Inputs format**:
    ```matlab
    SIZE = [6,7,8];
    
    % this ...
    [~,~,~,v1,v2,v3] = BH_multi_gridCoordinates(SIZE, 'Cartesian', 'cpu', {'none'}, 0, 1, 0);
    % ... is equivalent to this:
    [v1,v2,v3] = EMC_coordVectors(SIZE, 'cpu', {});  
    ```

2. **flgFreqSpace** -vs- **normalize**: flgFreqSpace triggers whether or not the vectors should be normalized by their number of elements, which is exactly what 'normalize' does.
    ```matlab
    % this ...
    [~,~,~,v1,v2,v3] = BH_multi_gridCoordinates(SIZE, 'Cartesian', 'cpu', {'none'}, ...
                                                1,... % flgFreqSpace
                                                1,... % flgShiftOrigin
                                                0);   % flgRad
    % ... is equivalent to this:
    [v1,v2,v3] = EMC_coordVectors(SIZE, 'cpu', {'normalize' true});  
    ```

3. **flgShiftOrigin** -vs- **origin**: As shown before, these control the origin convention.
    ```matlab
    % this ...
    [~,~,~,v1,v2,v3] = BH_multi_gridCoordinates(SIZE, 'Cartesian', 'cpu', {'none'}, ...
                                                0,... % flgFreqSpace
                                                1,... % flgShiftOrigin
                                                0);   % flgRad
    % ... is equivalent to this...
    [v1,v2,v3] = EMC_coordVectors(SIZE, 'cpu', {});  % or this:
    [v1,v2,v3] = EMC_coordVectors(SIZE, 'cpu', {'origin', 1});
    ```

4. **Shifts**: This is where we start to see that BH_multi_gridCoordinates was not designed to output *only* vectors.
    ```matlab
    % To apply a shift to the vectors with BH_multi_gridCoordinates,
    % the easiest is to use the 'single' transformation ('gridVector' does NOT work).
    SIZE = [6,7,8];
    shifts = [1;2;3];  % BH_multi_gridCoordinates needs a column vector
    
    % this ...
    [~,~,~,v1,v2,v3] = BH_multi_gridCoordinates(SIZE, 'Cartesian', 'cpu',
                                                {'single', [1,0,0;0,1,0;0,0,1], shifts, 'inv', 1, 1}, ...
                                                0,1,0);
    % ... is equivalent to this:
    [v1,v2,v3] = EMC_coordVectors(SIZE, 'cpu', {'shift', shifts});
    
    % note: EMC_coordVectors does not have directions ('inv'|'fwd') but the equivalent 
    %       is direction = 'inv' or 'forwardVector'.
    % note: EMC_coordVectors supports column and row vectors for SIZE and shifts.
    ```  

5. **halfGrid**:
    ```matlab
    % this ...
    [~,~,~,v1,v2,v3] = BH_multi_gridCoordinates(SIZE, 'Cartesian', 'cpu', {'none'}, 0,1,0, {'halfGrid'};
    % ... is equivalent to this:
    [v1,v2,v3] = EMC_coordVectors(SIZE, 'cpu', {'half', true});
    ```

6. **NaNs** and **empty** dimensions:
    ```matlab
    % Empty dimensions are not allowed.
    [v1,v2,v3] = EMC_coordVectors([0,10], 'cpu', {});  % will raise error with id=EMC:SIZE
    
    % unit dimensions (size=1) output a vector equals to NaN (easier to catch|identify than 1).
    [v1,v2,v3] = EMC_coordVectors([1,10,1], 'cpu', {});
    % v1 = v3 = NaN
    ```
    
#### EMC_coordGrids
1. Change the coordinate **SYSTEM**:
    ```matlab
    % this ...
    [g1,g2,g3,v1,v2,v3] = BH_multi_gridCoordinates(SIZE, 'Cartesian', 'cpu', {'none'}, 0,1,0);
    % ... is equivalent to this:
    [g1,g2,g3,v1,v2,v3] = EMC_coordGrids('Cartesian', SIZE, 'cpu', {});
    
    % note: EMC_coordGrids uses strcmpi to check SYSTEM, so 'cartesian works too'.
    
    % note: EMC_coordGrids has extacly the same OPTION parameters than EMC_coordVectors
    help EMC_coordGrids
    help EMC_coordVectors
    ```

2. **Radial** grid:
    ```matlab
    % this ...
    g1 = BH_multi_gridCoordinates(SIZE, 'Cartesian', 'cpu', {'none'}, ...
                                  0,... % flgFreqSpace
                                  1,... % flgShiftOrigin
                                  1);   % flgRad
    % ... is equivalent to this:
    g1 = EMC_coordGrids('Radial', SIZE, 'cpu', {});
    % note: techically radial grids are cartesian grids, but for simplicity,
    %       SYSTEM accepts 'Radial' in addition of 'cartiesian', 'cylindrical' and 'spherical'.
    
    % To compute radial grids with a given radius, I [TF] use EMC_coordVectors:
    radius = [30,20];
    SIZE = [100,100];
    [v1,v2] = EMC_coordVectors(SIZE, 'cpu', {});
    radialGrid = (v1'./radius(1)).^2 + (v2./radius(2)).^2;
    % note: v1 is transposed to broadcast the vectors into a grid.
    ```

#### EMC_coordTransform
TODO: add this part, but not a priority, since the transformation are now done with mex files.
