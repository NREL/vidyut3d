python get_error_1d.py case10 0.1 1 > errconv5
python get_error_1d.py case20 0.1 1 >> errconv5
python get_error_1d.py case40 0.1 1 >> errconv5
python get_error_1d.py case80 0.1 1 >> errconv5
python get_error_1d.py case160 0.1 1 >> errconv5
gnuplot pic_adv_weno.gp
