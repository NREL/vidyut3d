./*.ex inputs vidyut.dt=0.01 vidyut.num_timestep_correctors=4
python get_L2norm_error.py plt00001 0 0.1 0.01 > err4
./*.ex inputs vidyut.dt=0.005 vidyut.num_timestep_correctors=4
python get_L2norm_error.py plt00001 0 0.1 0.005 >> err4
./*.ex inputs vidyut.dt=0.0025 vidyut.num_timestep_correctors=4
python get_L2norm_error.py plt00001 0 0.1 0.0025 >> err4
./*.ex inputs vidyut.dt=0.00125 vidyut.num_timestep_correctors=4
python get_L2norm_error.py plt00001 0 0.1 0.00125 >> err4
