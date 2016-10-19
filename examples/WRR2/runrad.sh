#@echo off
#rem normally run a whole year, for this example only 10 days!
python ../../e2o_dstools/e2o_radiation.py -S 1 -E 10 -D gtopo05min.map -d WRR2 -O outputrad -s 5 -e 22 -l DEBUG
