# Synt_v3

Modified pipeline for DNAWorks tool.

First of all you need install xcode for macOS 

`xcode-select --install`

After that clone repo with DNAWorks

`git clone https://github.com/davidhoover/DNAWorks.git`

Go to the directory with DNAWorks

`cd DNAWorks`
`make`

If you receve an error when code is compiling use command 

`export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"`

and try again 
`make`

Go back to folder

`cd ..`

Install virtual environment with comda/micromamba

`conda env create -f env.yml`
`conda activate Synt_v3`

`python init.py -n My_Test_Project MAQLGKLLKEQKYDRQLRLWGDHGQEALESAHVCLINATATGTEILKNLVLPGIGSFTIIDGNQVSGEDAGNNFFLQRSSIGKNRAEAAMEFLQELNSDVSGSFVEESPENLLDNDPSFFCRFTVVVATQLPESTSLRLAD`

You'll get new folder with your project name.
Folder contains files with generated sequences, primers lists, plots with heterodimer Tm, GC-content and so on.
