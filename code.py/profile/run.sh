# by code line (inaccurate)
for py in model/*.py; do sed -i 's/\#\@profile/\@profile/' $py; done
(date; kernprof -lv -u 1 profile/run.py) > profile/profile.l.txt
for py in model/*.py; do sed -i 's/\@profile/\#\@profile/' $py; done
rm run.py.lprof
# by code element (accurate)
(date; kernprof -v -u 1 profile/run.py) > profile/profile.e.txt
rm run.py.prof
