# Rabies Minion Amplicon Analysis Pipeline 

This pipeline takes basecalled, demultiplexed, and trimmed fastq files from Nanopore sequencing.
It QC the fastq files using fastqc and aggregates the reports using multiqc.
It mapps to the reference rabies genomes to get coverage of N and G area.
It uses Canu to correct the reads into amplicons.
The amplicons are clustered using CD-hit-est. 
Clusters with greater than N amplicons are selected, and
then, polished by Medaka the fastq files again.

This workflow depends on Bioconda to install the environmental dependencies.
And then run the attached script using the environment.

## Summary - Installation 
1. Clone Repository 
2. Install Conda if not already in environment
3. Create conda environment

## Summary - How to run after installation.
1. Activate conda environment - `conda activate RabiesMAAP`
2. Run the GUI - `~/RabiesMAAP/RabiesMAAP_GUI.py`
4. Test data in test folder

## Clone this code using GIT

### Install git for Debian systems using the following command (if necessary)
```
sudo apt update
sudo apt install git
```

##Installation directions 
These instructions install the code into your home path. Change the instructions if appropriate. 

### Clone the code from repository
```
cd ~
git clone https://github.com/jiangweiyao/RabiesMAAP.git
```

### Install Both Miniconda and Environment (Skip next 2 steps)
You can run the prepackaged script install_all.sh to install both Miniconda and make the Conda environment into your home directory (recommended) by using the following command
```
. ~/RabiesMAAP/install_miniconda.sh
```


### Install Miniconda (if no Conda is install on system). 
You can run the prepackaged script install_miniconda.sh to install into your home directory (recommended) by using the following command
```
. ~/RabiesMAAP/install_miniconda.sh
```

Detailed instruction on the the Miniconda website if anything goes wrong:
https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

### Clone the environment. Need to do once.

We use conda to create an environment (that we can activate and deactivate) to install our dependent software and resolve their dependencies. This environment is called "RabiesMAAP". The following command assumes your environment file is in your home path. Modify as appropriate.

```
conda env create -f ~/RabiesMAAP/environment.yml
```

The command to generate the environment originally is in the included `conda_env.txt` file. 

## Run the code.

Activating your environment makes the software you installed in that environment available for use. You will see "(RabiesMAAP)" in front bash after activation.
```
conda activate RabiesMAAP
```

When you are finished running the workflow, exit out of your environment by running `conda deactivate`. Deactivating your environment exits out of your current environment and protects it from being modified by other programs. You can build as many environments as you want and enter and exit out of them. Each environment is separate from each other to prevent version or dependency clashes. The author recommands using Conda/Bioconda to manage your dependencies.

## Known issues
### Issue 1: Perl clashing for canu
Some systems preset the PERL5LIB variable, causing canu to look for the perl5 library in that location rather than the default location in the conda environment. You can manually reset the PERL5LIB variance using the following command
`PERL5LIB=`

## Author
Jiangwei Yao
Crystal Gigante

## License 
Apache License 2.0
