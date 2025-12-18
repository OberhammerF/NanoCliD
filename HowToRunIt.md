# How to run the NanoCliD pipeline

1. suffer
2. start interactive job 
    `srun -p gpu --cpus-per-task=8 --mem=16G --gpus-per-node=1 --time=8:00:00 --pty bash -`
3. load miniconda 24.9.2?
    `module load miniconda/24.9.2`
4. activate nanoclid environment
    `conda activate nanoclid`
5. activate venv 
    1. `cd /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD`
    2. `source venv/bin/activate`
6. set the correct tmp dir
    1. `export TMPDIR=/tmp`

7. run pipeline
    `python nanoclid_FO.py test   -p standalone   -R /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/nanoclid_ref/  -C /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD/singularity -M /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD/annotations/dorado_models/ -X guppy  &> pipeline.log`

    you need a hg19.fa, the index, a dictionary (made with picardtools), and a genome file. maybe even more