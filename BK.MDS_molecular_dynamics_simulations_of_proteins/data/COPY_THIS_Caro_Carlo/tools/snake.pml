reset

show sticks
set stick_radius, 1.5
alter all, vdw=1.25
bg_color white


select aa, index 1+7+13
alter aa, vdw=3.0
rebuild
show spheres, aa

set movie_fps, 60

spectrum
deselect

zoom center, 45

set specular_intensity, 0.0
unset depth_cue
