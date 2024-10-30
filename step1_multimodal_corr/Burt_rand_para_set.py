from brainsmash.workbench.geo import volume
from brainsmash.mapgen.eval import sampled_fit
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio

def volume_gen(coord_file,output_dir,brain_map,ns,knn,output_pic_name):
    filenames = volume(coord_file,output_dir)
    kwargs = {'ns': ns,
          'knn': knn,
          'pv': 70,
          }
    
    fig, rho = sampled_fit(brain_map, filenames['D'], filenames['index'], nsurr=1000, **kwargs)
    fig.savefig(output_dir + '\\' + output_pic_name + '_' + str(rho) + '.png')
    
    return

