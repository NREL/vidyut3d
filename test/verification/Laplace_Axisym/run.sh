mpirun -n 2 ./*.ex inputs2d max_step=2 vidyut.pot_bc_lo = 1 2 vidyut.pot_bc_hi = 1 2 prob.left_dirc=1 prob.right_dirc=1 prob.left_neumann=0 prob.right_neumann=0 
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh

mpirun -n 2 ./*.ex inputs2d max_step=2 vidyut.pot_bc_lo = 4 2 vidyut.pot_bc_hi = 1 2 prob.left_dirc=1 prob.right_dirc=1 prob.left_neumann=0 prob.right_neumann=0 
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh

mpirun -n 2 ./*.ex inputs2d max_step=2 vidyut.pot_bc_lo = 4 2 vidyut.pot_bc_hi = 1 2 prob.left_dirc=0 prob.right_dirc=1 prob.left_neumann=1 prob.right_neumann=0 
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh

mpirun -n 2 ./*.ex inputs2d max_step=2 vidyut.pot_bc_lo = 3 2 vidyut.pot_bc_hi = 1 2 prob.left_dirc=0 prob.right_dirc=1 prob.left_neumann=1 prob.right_neumann=0 
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh

mpirun -n 2 ./*.ex inputs2d max_step=2 vidyut.pot_bc_lo = 1 2 vidyut.pot_bc_hi = 4 2 prob.left_dirc=1 prob.right_dirc=1 prob.left_neumann=0 prob.right_neumann=0 
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh

mpirun -n 2 ./*.ex inputs2d max_step=2 vidyut.pot_bc_lo = 1 2 vidyut.pot_bc_hi = 4 2 prob.left_dirc=1 prob.right_dirc=0 prob.left_neumann=0 prob.right_neumann=1 
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh

mpirun -n 2 ./*.ex inputs2d max_step=2 vidyut.pot_bc_lo = 1 2 vidyut.pot_bc_hi = 3 2 prob.left_dirc=1 prob.right_dirc=0 prob.left_neumann=0 prob.right_neumann=1
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh

mpirun -n 2 ./*.ex inputs2d amr.max_level=1 max_step=2 vidyut.pot_bc_lo = 1 2 vidyut.pot_bc_hi = 1 2 prob.left_dirc=1 prob.right_dirc=1 prob.left_neumann=0 prob.right_neumann=0
~/fextract -s pot1d.dat -d 0 -v Potential plt00001
python get_error.py plt00001 pot1d.dat
. ../../clean.sh
