## BH_bandpass3d &rarr; EMC_getBandpass.m
These functions are almost indentical.

### Signature:
```matlab
BANDPASS = BH_bandpass3d(SIZE, HIGH_THRESH, HIGHPASS, LOWPASS, METHOD, PIXEL, varargin);
% vs
BANDPASS = EMC_getBandpass(SIZE, PIXEL, HIGHPASS, LOWPASS, METHOD, OPTION);
```

### Noticeable difference:
1. **Gaussian roll**: This is the main difference between ```BH_bandpass3d``` and ```EMC_getBandpass``` at the moment. The default with ```EMC_getBandpass``` consists of an adjustable roll off whereas ```BH_bandpass3d``` has a fixed ~7 pixels size. The default of ```EMC_getBandpass``` will soon be changed to the same fixed roll off as ```BH_bandpass3d```.
  ```matlab
  % this ...
  bandpass = BH_bandpass3d([100,100], 0, 30, 14, 'cpu', 1);
  % ... is 'equivalent' to this:
  bandpass = EMC_getBandpass([100,100], 1, 30, 14, 'cpu', {});
  
  % note: for a lowpass filter -> HIGHPASS should be NaN or 0.
  % note: for a highpass filter -> LOWPASS should be Nan or 0.
  ```

2. **SIRT-like** lowpass: With ```BH_bandpass3d```, if ```LOWPASS < 0```, the gaussian roll is from ```abs(LOWPASS)``` to Nyquist (threshold at Nyquist = 1e-3). These types of extended roll off are fully supported in ```EMC_getBandpass```; see example below.
  ```matlab
  % this ...
  bandpass = BH_bandpass3d([100,100], 0, 30, -14, 'cpu', 1);
  % ... is 'equivalent' to this:
  bandpass = EMC_getBandpass([100,100], 1, 30, 14, 'cpu', {'lowpasRoll', 'extended'});
  
  % To see all the optional parameters related to the gaussian roll:
  help EMC_getBandpass
  % One can modify both the gaussian std and the threshold, for the lowpass cut and highpass cut.
  
  % note: if the extended roll off is shorter than the default roll off because
  %       LOWPASS is too close to Nyquist, the function switch back to the default taper.
  %       This is also valid for HIGHPASS and the DC component.
  ```

3. **nyquistHigh**: If ```PIXEL``` is not numeric, ```BH_bandpass3d``` computes a lowpass filter with the ```LOWPASS``` at nyquist. If ```PIXEL == 'nyquistHigh'```, the ```highCut``` is adjusted to allow for at least a 7 pixel taper. ```highCut``` is also adjusted in normal cases to allow for at least a 7 pixel taper at the zero frequency. This behavior is not present in ```EMC_getBandpass``` yet.
  ```matlab
  % lowpass filter at Nyquist
  bandpass = EMC_getBandpass([100,100], 2, nan, 4, 'cpu', {});  % or
  bandpass = EMC_getBandpass([100,100], 2, nan, 'nyquist', 'cpu', {});
  ```

4. **half (non-redundant) grids**:
  ```matlab
  % this ...
  bandpass = BH_bandpass3d([100,100], 0, 30, 14, 'cpu', 1, 'halfgrid');
  % ... is 'equivalent' to this:
  bandpass = EMC_getBandpass([100,100], 1, 30, 14, 'cpu', {'half', true});
  ```

5. **Centered grids**: By default both ```BH_bandpass3d``` and ```EMC_getBandpass``` compute not-centered grids (```origin=-1```). ```EMC_getBandpass``` has access to every supported ```'origin'``` and can therefore compute centered grid directly.
  ```matlab
  bandpass = EMC_getBandpass([100,100], 1, 30, 14, 'cpu', {'origin', 1});  % origin = -1, 0, 1, or 2
  % of course, centered half grids are also possible:
  bandpass = EMC_getBandpass([100,100], 1, 30, 14, 'cpu', {'origin', 1; 'half', true});
  ```

