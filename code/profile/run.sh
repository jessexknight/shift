# WARNING: we modify original files: don't tinker with the 's/find/replace/'!
for py in model/*.py; do sed -i 's/\#\@profile/\@profile/' $py; done
(date; kernprof -lv -u 1 profile/run.py) > profile/profile.txt
for py in model/*.py; do sed -i 's/\@profile/\#\@profile/' $py; done
rm run.py.lprof
