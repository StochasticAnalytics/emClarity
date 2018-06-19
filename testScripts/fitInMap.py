# Script to fit the two half-sets together using chimera's fit in map.
#

# Only a local optimization is done so the initial position must be close
# to the correct fit.
#
# map arg1 is fit into map arg2
# arg3 is a string prefix to name the xform file output
#
# Use the inverse of the rotMat and the translation as given with 
# bah forward convention. The origin must be set in the headers to be centered, 
# which seems to be floor(nX/2) 

from sys import argv

def fit_map_in_map(map1_path, map2_path, xformName,
                   initial_map1_transform = None,
                   map1_threshold = None,
                   ijk_step_size_min = 0.01,    # Grid index units
                   ijk_step_size_max = 0.5,     # Grid index units
                   max_steps = 5000,
                   optimize_translation = True,
                   optimize_rotation = True):


  # Files have to have file suffix indicating volume format.
  from VolumeViewer import open_volume_file
  map1 = open_volume_file(map1_path)[0]  # Assume files contain just one array
  map2 = open_volume_file(map2_path)[0]

  if initial_map1_transform:
    from Matrix import chimera_xform
    xf = chimera_xform(initial_map1_transform)
    map1.surface_model().openState.globalXform(xf)
    
  use_threshold = (map1_threshold != None)
  
  from FitMap.fitmap import map_points_and_weights, motion_to_maximum
  points, point_weights = map_points_and_weights(map1, use_threshold)

  if len(points) == 0:
    if use_threshold:
      print 'No grid points above map threshold.'
    else:
      print 'Map has no non-zero values.'
    return

  move_tf, stats = motion_to_maximum(points, point_weights, map2, max_steps,
                                     ijk_step_size_min, ijk_step_size_max,
                                     optimize_translation, optimize_rotation)

  import Matrix
  if initial_map1_transform:
    move_tf = Matrix.multiply_matrices(move_tf, initial_map1_transform)

  header = ('\nFit map %s in map %s using %d points\n'
            % (map1.name, map2.name, stats['points']) +
            '  correlation = %.4g, overlap = %.4g\n'
            % (stats['correlation'], stats['overlap']) +
            '  steps = %d, shift = %.3g, angle = %.3g degrees\n'
            % (stats['steps'], stats['shift'], stats['angle']))
  print header
  f = open(xformName,'w')
  for i in range(4):  
    for j in range(3):
     # print move_tf[j][i]
      f.write('{:f}'.format(move_tf[j][i]))
      f.write("\n")

  f.close()
  print move_tf
  tfs = Matrix.transformation_description(move_tf)
  print tfs


# -----------------------------------------------------------------------------
# Example call fitting map into itself.
#

#p1 = 'rotbe.mrc'
#p2 = 'rotae.mrc'
t = fit_map_in_map(argv[1],argv[2],argv[3])


