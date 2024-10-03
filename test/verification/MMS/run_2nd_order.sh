./*.ex inputs_x amr.max_level=0 amr.n_cell=8 4 4 vidyut.hyp_order=2 vidyut.dt=0.001 max_step=100000
python get_L2norm_error.py plt00001 > err2
./*.ex inputs_x amr.max_level=0 amr.n_cell=16 4 4 vidyut.hyp_order=2 vidyut.dt=0.001 max_step=100000
python get_L2norm_error.py plt00001 >> err2
./*.ex inputs_x amr.max_level=0 amr.n_cell=32 4 4 vidyut.hyp_order=2 vidyut.dt=0.001 max_step=100000
python get_L2norm_error.py plt00001 >> err2
./*.ex inputs_x amr.max_level=0 amr.n_cell=64 4 4 vidyut.hyp_order=2 vidyut.dt=0.001 max_step=100000
python get_L2norm_error.py plt00001 >> err2
./*.ex inputs_x amr.max_level=0 amr.n_cell=128 4 4 vidyut.hyp_order=2 vidyut.dt=0.001 max_step=100000
python get_L2norm_error.py plt00001 >> err2
