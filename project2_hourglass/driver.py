from Animation import Container_Animation
from Initialize import ParticleInitialize
import matplotlib.pyplot as plt
########## PARAMS ###############
num_frames = 20
delay      = 0
vid_format = ".mp4"

run_backwards  = False
save_animation = False
print_frame    = True
tile_domain    = False
crunch_domain  = False
animate_it     = True
show_figs      = False

initializer = ParticleInitialize()
c,distance_matrix,force,integrate,xlim,ylim,case = initializer()

print c

if animate_it:
  animator = Container_Animation(c,integrate,force,xlim,ylim,save_animation=save_animation,print_frame=print_frame,filename=case,num_frames=num_frames)
  if save_animation:
    animator.save()
  else:
    animator.show()

# for i in range(num_frames):
#   if not( i%100 ):
#     print i
#   c.integrate()

#from datetime import datetime
#now = datetime.now()
#datestr = now.strftime("%Y-%m-%d_%H:%M") 
#
#import pickle
#f = open(datestr+"force_run.dump","w")
#pickle.dump([max_forces,run_dat],f)
#f.close()
#plt.show()
