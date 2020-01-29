## BH_padZeros3d &rarr; EMC_resize
### Signature:
```matlab
 [OUT_IMG] = EMC_resize(IN_IMG, LIMITS, OPTION);
 [OUT_IMG] = BH_padZeros3d(IN_IMG, PADLOW, PADTOP, METHOD, PRECISION, varargin);
```

### Not supported:
1.  PADLOW can be 'fwd' or 'inv' for some reason. This is not supported in EMC_resize.

### Noticeable differences:
1. **OPTION**:
  As every EMC functions, **EMC_resize** uses the **OPTION** parameter to deal with parameters with default values. **OPTION** can be a cell or a structure:
	 ```matlab
	 OPTION1 = {param1, value1; param2, value2};
	 % is equivalent to
	 OPTION2.param1 = value1;
	 OPTION2.param2 = value2;
     ```
     
2. **METHOD**:
  **EMC_resize** does NOT allow the input image to change 'METHOD' (switch cpu <-> gpu). If one wants to change the method, one needs to explicitly push to device or gather to host before calling **EMC_resize**. 
	 ```matlab
	 % This is equivalent - push to device
	 img = ones(50,50);
	 [imgGpu] = EMC_resize(gpuArray(img), [0,0,0,0], {});
	 [imgGpu] = BH_padZeros3d(img, [0,0], [0,0], 'GPU', 'double');
	 
	 % This is equivalent - gather to host
	 img = ones(50,50, 'gpuArray');
	 [imgCpu] = EMC_resize(gather(img), [0,0,0,0], {});
	 [imgCpu] = BH_padZeros3d(img, [0,0], [0,0], 'cpu', 'double');
	 
	 % EMC_resize does NOT have a METHOD parameter. The METHOD will be deduced
	 % from the input IMAGE. This is to prevent non-explicit transfer between
	 % host and device, as in BH_padZeros3d. Moreover, the METHOD is not necessary
	 % respected in BH_padZeros3d (ex: if input is 'cpu' and METHOD='GPU', the taper
	 % is applied on the input, so on the cpu).
   ```

3. **PADLOW** and **PADTOP** are merged into **LIMITS**:
	 ```matlab
	 padlow = [1,2];  % [xlow, ylow]
	 padtop = [3,4];  % [xhigh, yhigh]
	 % is equivalent to
	 limits = [1,3,2,4]  % [xlow, xhigh, ylow, yhigh]
	 
	% Example - This equivalent
	 [out1] = EMC_resize(img, limits, {});  % by default, if padding, apply taper. See below for more details.
	 [out2] = BH_padZeros3d(img, padlow, padtop, 'cpu', 'singleTaper')
   ```

4. **Tapers**:
  **EMC_resize** tapers include the padding value in their last pixel (they're one pixel longer by default).
	 ```matlab
	 % |---------img(ones)-taper(cosine)-------| pad(zeros)
	   [...,1,  0.95,0.81,0.61,0.39,0.19,0.05,   0,0,... ]  % BH_padZeros3d
	 
	 % |---------img(ones)-taper(cosine)---------| pad(zeros)
	   [...,1,  0.95,0.81,0.61,0.39,0.19,0.05,0,   0,0,... ]  % EMC_resize
	   % note the extra 0 at the end of taper.
     ```
	To generate the same roll off, **EMC_resize** will use an additional pixel of the input image. In counterpart, it makes things more predictable and intuitive when the image is not padded (see 5.example2). Moreover, **EMC_resize** allow the user to use its own taper:
	 ```matlab
	 % change the size of the taper to 20 pixels.
	 [img2dDouble] = EMC_resize(img, [0,0,0,0], {'taper', {'cosine', 20}});
	 % switch to linear taper
	 [img2dDouble] = EMC_resize(img, [0,0,0,0], {'taper', {'linear', 20}});
	 % use your own taper
	 taper = [1, 0.8, 0.6, 0.4];
	 [img2dDouble] = EMC_resize(img, [0,0,0,0], {'taper', taper);
	 
	 % EMC_resize is calling EMC_taper to create tapers.
	 help EMC_taper
   ```

5. **PRECISION**:
  In **BH_padZeros3d**, **PRECISION** manages the precision of the output and if a taper is meant to be applied to the input image. In  **EMC_resize**, this responsibility is shared between 3 optional parameters ('**precision**', '**taper**' and '**force_taper**'). By default, the precision of the output image is the same as the input image. If you want to change it, add the following to OPTION: ```{'precision', desiredPrecision}```. To turn off the taper: ```{'taper', false}```. To tape the image on every edges whether or not they are padded|cropped: ```{'force_taper', true}```.

   #####  EXAMPLE 1: no padding, no cropping, no taper, but cast from single to double:
   ```matlab
	imgSingle = ones(50,50,'single');
	
	% This is equivalent:
	[imgDouble1] = EMC_resize(imgSingle, [0,0,0,0], {'precision', 'double'});
	[imgDouble2] = BH_padZeros3d(imgSingle, [0,0], [0,0], 'cpu', 'double');
   ```

   #####  EXAMPLE 2: Apply a taper to every edge of the image, without cropping or padding:
   ```matlab
	% This is equivalent - apply a taper that goes to zeros (default)
	[imgSingle1] = EMC_resize(imgSingle, [0,0,0,0], {'force_taper', true});
	[imgSingle2] = BH_padZeros3d(imgSingle, [0,0], [0,0], 'cpu', 'singleTaper');
	 
	% This is equivalent - apply a taper that goes to 5
	[imgSingle1] = EMC_resize(imgSingle, [0,0,0,0], {'force_taper', true});
	                                                 'value', 5});
	[imgSingle2] = BH_padZeros3d(imgSingle, [0,0], [0,0], 'cpu', 'singleTaper', 5);
	
	% In this case (no padding), BH_padZeros3d outputs can be surprising because the 
	% taper doesn't go to the desired value 5. It makes more sense when the image is
	% padded though.
	% As shown previously, EMC_resize default tapers 'include the last pixel', which
	% makes the edges equal to the padding value (with or without padding|cropping).
   ```

   #####  EXAMPLE 3: default padding and taper convention
   ```matlab
	% By default, EMC_resize will tape any edges that are padded. To deactivate
	% the taper, use:
	[imgSingle1] = EMC_resize(imgSingle, [10,20,30,40], {'taper', false});
	[imgSingle2] = BH_padZeros3d(imgSingle, [10,30], [20,40], 'cpu', 'single');
	
	% 
	% This is equivalent - apply taper with padding...
	[imgSingle1] = EMC_resize(imgSingle, [1,2,3,4], {});
	[imgSingle2] = BH_padZeros3d(imgSingle, [1,3], [2,4], 'cpu', 'singleTaper')
	 
	 % ... but this is not, because EMC_resize applies (by default) a taper only to the
	 % edges that are padded (here xlow is not padded, so no taper on this side).
	 [imgSingle1] = EMC_resize(imgSingle, [0,2,3,4], {});
	 [imgSingle2] = BH_padZeros3d(imgSingle, [0,3], [2,4], 'cpu', 'singleTaper')
	 
	 % This is where 'force_taper' is useful: if true, it applies the taper to every
	 % edges, whether or not they are padded|cropped. This is now equivalent:
	 [imgSingle1] = EMC_resize(imgSingle, [0,2,3,4], {'force_taper', true});
	 [imgSingle2] = BH_padZeros3d(imgSingle, [0,3], [2,4], 'cpu', 'singleTaper')
    ```

6. **Padding value**:
	 ```matlab
	 % This is equivalent - pad with value = 12
	 [out1] = EMC_resize(img, [10,20,30,40], {'value', 12});
	 [out2] = BH_padZeros3d(img, [10,30], [20,40], 'cpu', 'singleTaper', 12);
	 
	 % This is equivalent - pad with white gaussian noise (scaled on img mean and std)
	 [out1] = EMC_resize(img, [10,20,30,40], {'value', 'uniform'});
	 [out2] = BH_padZeros3d(img, [10,30], [20,40], 'cpu', 'singleTaper', 'uniform');
	 % technically, BH_padZeros3d will treat everything that is not numeric
	 % the same way, so '' will work too...
	 
	 % Btw, to pad with the mean, the value parameter accepts 'mean', so
	 [out1] = EMC_resize(img, [10,20,30,40], {'value', mean(img(:))});
	 % is equivalent to
	 [out1] = EMC_resize(img, [10,20,30,40], {'value', 'mean'});
   ```


7. **fourierOverSample** vs '**origin**':
  Every EMC functions use the same parameter to deal with the origin|center of an array (1d,2d or 3d): optional parameter 'origin'. origin = -1 means the input array is meant to be treated as non-centered (zero frequency first). It is equivalent to fourierOverSample=1.
	 ```matlab
	 % This is equivalent
	 [img2dDouble] = EMC_resize(img, [10,20,30,40], {'origin', -1});
	 [img2dDouble] = BH_padZeros3d(img, [10,30], [20,40], 'cpu', 'singleTaper', 1,1);
	 % I am not sure how you [Ben] originaly planned this, but if 
	 % len(varargin) > 1, fourierOverSample=1, so varargin = {1,1} works.
   ```
