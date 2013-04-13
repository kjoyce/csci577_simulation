from Animation import Container_Animation
from Initialize import ParticleInitialize
import matplotlib.pyplot as plt
########## PARAMS ###############
num_frames = 8000
ekg_length = 80
delay      = 0
vid_format = ".mp4"

run_backwards  = False
save_animation = False
print_frame    = False
tile_domain    = False
crunch_domain  = False
animate_it     = False
show_figs      = False

max_forces = dict()
run_dat = []
for n_sled in (9,13,17):
  for load in range(-20,45,5):
    if n_sled == 9 and load < 0:
      continue
    initializer = ParticleInitialize(n_sled,load)
    c,distance_matrix,force,integrate,xlim,ylim,pull_force_lim,ave_vel_lim,case = initializer()

    case += '_time{}'.format(num_frames)

    filename = case
    if animate_it:
      animator = Container_Animation(c,integrate,force,xlim,ylim,pull_force_lim,ave_vel_lim,run_backwards=run_backwards,save_animation=save_animation,print_frame=print_frame,tile_domain=tile_domain,crunch_domain=crunch_domain,vid_format=vid_format,filename=filename,num_frames=num_frames,ekg_length=ekg_length)
      if save_animation:
	animator.save()
      else:
	animator.show()

    for i in range(num_frames):
      if not( i%100 ):
	print i
      c.integrate()

    fmax = max(c.pull_force)
    max_forces[case] = fmax
    plt.figure()
    plt.plot(range(num_frames),c.pull_force[0:num_frames])
    plt.title("Pull Force vs. Time")
    plt.xlabel("Time (.01) (non-dimensionalized)")
    plt.ylabel("Pull Force (non-dimensionalized)")
    plt.xlim((-15,45))
    plt.savefig('pullforce_'+case+'_eqscale.pdf')
#    plt.figure()
#    plt.plot(range(num_frames),c.avg_velocities[0:num_frames])
#    plt.title("Average Horizontal Velocity of Structure vs. Time")
#    plt.xlabel("Time (.01) (non-dimensionalized)")
#    plt.ylabel("Ave. x-velocity (non-dimensionalized)")
#    plt.savefig('avevelocity_'+case+'.pdf')
    if show_figs:
      plt.show()
    run_dat.append([c,distance_matrix,force,integrate,xlim,ylim,pull_force_lim,ave_vel_lim,case]) 

from datetime import datetime
now = datetime.now()
datestr = now.strftime("%Y-%m-%d_%H:%M") 

import pickle
f = open(datestr+"force_run.dump","w")
pickle.dump([max_forces,run_dat],f)
f.close()
#plt.show()
