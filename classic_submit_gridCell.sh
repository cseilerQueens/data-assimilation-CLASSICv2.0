#!/bin/sh
#SBATCH --account=def-cseiler-ab 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# #SBATCH --exclusive
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --job-name=test01
#SBATCH --output=/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/classic.out
#SBATCH --error=errors_CLASSIC
# #SBATCH --mail-user=christian.seiler@queensu.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

cd /home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/test01
/home/cseiler/CLASSICv2.0/classic/bin/CLASSIC_parallel_intel job_options.txt 250/60
touch /home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/classicFinished
#!/bin/sh
#SBATCH --account=def-cseiler-ab 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# #SBATCH --exclusive
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --job-name=test01
#SBATCH --output=/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/classic.out
#SBATCH --error=errors_CLASSIC
# #SBATCH --mail-user=christian.seiler@queensu.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

cd /home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/test01
/home/cseiler/CLASSICv2.0/classic/bin/CLASSIC_parallel_intel job_options.txt 250/60
touch /home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation-CLASSICv2.0/simulations/classicFinished
