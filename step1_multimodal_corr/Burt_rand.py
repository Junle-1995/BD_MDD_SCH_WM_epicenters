from brainsmash.workbench.geo import volume
from brainsmash.mapgen.eval import sampled_fit
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio

def volume_gen(coord_file,output_dir,brain_map,ns,knn,output_file_name):
    filenames = volume(coord_file,output_dir)
    kwargs = {'ns': ns,
          'knn': knn,
          'pv': 70,
          }
    
    gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], **kwargs)
    surrogate_maps = gen(10000)

    sio.savemat(output_dir + '\\' + output_file_name,{'surrogate_maps':surrogate_maps})
    return

