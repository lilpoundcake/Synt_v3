# Synt_v3

Modified pipeline for DNAWorks tool.

## MacOS

First of all you need install xcode for macOS 

`xcode-select --install`

After that clone repo with DNAWorks

`git clone https://github.com/davidhoover/DNAWorks.git`

Go to the directory with DNAWorks

`cd DNAWorks`

`make`

If you receve an error when code is compiling use command 

`export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"`

and try again:

`make`

Go back to folder

`cd ..`

## Linux

Install Gfortran:

`apt-get install gfortran`

## Conda environment

If you steel don't have [Miniconda](https://docs.anaconda.com/miniconda/) or [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) - The time has come. Go to site and find appropriate information.

Install virtual environment with conda/micromamba

`conda env create -f env.yml`

`conda activate Synt_v3`

For showing help enter:

`python init.py -h`

Tool supports several parameters:

- **-n** - name of project. Default name is *synt_gene*
- **-p** - *protein* or *nucleotide*. Default sequence type is *protein*
- **-e** - expression system. You may select one from list: *cho*, *ecoli*, *sf9*, *human*. Default is *cho*
- **-t** - preferable anneling temperature for PCR. Default temperature is 65 degrees Celsius.
- **-l** - by default primers will be designed with length less than 50 nt. You can choose another. Optimal value is between 50 and 70 nt.

---

Let's try to design oligos:

`python init.py -n My_Test_Project -p protein -e cho MAQLGKLLKEQKYDRQLRLWGDHGQEALESAHVCLINATATGTEILKNLVLPGIGSFTIIDGNQVSGEDAGNNFFLQRSSIGKNRAEAAMEFLQELNSDVSGSFVEESPENLLDNDPSFFCRFTVVVATQLPESTSLRLAD`

You'll get new folder with your project name.

Folder contains files with generated sequences, primers lists, plots with heterodimer Tm, GC-content and so on.
