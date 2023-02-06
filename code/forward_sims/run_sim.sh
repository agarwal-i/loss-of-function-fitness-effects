# Get simulators from Fuller et al 2019

# compile simulators
## standard Shiffels Durbin for EUR
# g++ -O3 -std=c++11 BRand.cpp population.cpp main_sd_mod_allele_age.cpp -o simulator_mod_sd

# run simulations
for i in hs10_1; do qsub -l mem=2G,time=1:00:00 -cwd batch_sim.sh 100000 1e-6 0.01 0 0.2 $i test_fixedu sel; done
for i in hs10_2; do qsub -l mem=2G,time=1:00:00 -cwd batch_sim.sh 100000 1e-6 0.01 0 0.2 $i test_fixedu sel; done
for i in hs50_1; do qsub -l mem=2G,time=1:00:00 -cwd batch_sim.sh 100000 1e-6 0.01 0 1 $i test_fixedu sel; done
for i in hs50_2; do qsub -l mem=2G,time=1:00:00 -cwd batch_sim.sh 100000 1e-6 0.01 0 1 $i test_fixedu sel; done
for i in hs50_3; do qsub -l mem=2G,time=1:00:00 -cwd batch_sim.sh 100000 1e-6 0.01 0 1 $i test_fixedu sel; done
for i in hs50_4; do qsub -l mem=2G,time=1:00:00 -cwd batch_sim.sh 100000 1e-6 0.01 0 1 $i test_fixedu sel; done

##output
# cat out/XX/*txt > out/XX.txt | gzip
