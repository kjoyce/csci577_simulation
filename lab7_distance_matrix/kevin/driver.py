from Animation import Container_Animation
from Initialize import ParticleInitialize
########## PARAMS ###############
num_frames = 1000
ekg_length = 80
delay      = 0
vid_format = ".mp4"

run_backwards  = False
save_animation = False
print_frame    = False
tile_domain    = False
crunch_domain  = False

#case = "one"
#case = "two"
#case = "three"
#case = "four"
#case = "six"
#case = "eight"
#case = "line"
case = "square_lattice"
#case = "crunch_square_lattice"
#case = "triangle_lattice"
#case = "crunch_triangle_lattice"
#################################

initializer = ParticleInitialize(case)
c,distance_matrix,force,integrate,xlim,ylim,pot_energy_lim,kin_energy_lim,tot_energy_lim,pressure_lim = initializer()

filename = case
if crunch_domain:
  filename = "crunch_"+filename
animator = Container_Animation(c,integrate,force,xlim,ylim,pot_energy_lim,kin_energy_lim,tot_energy_lim,pressure_lim,run_backwards=run_backwards,save_animation=save_animation,print_frame=print_frame,tile_domain=tile_domain,crunch_domain=crunch_domain,vid_format=vid_format,filename=filename,num_frames=num_frames,ekg_length=ekg_length)
animator.show()
