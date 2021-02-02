# WORK IN PROGRESS
There are some bugs known to me, some aren't. Before getting desperate, please check out the Issues that are already opened and discussed. If you can't find your problem there, don't hesitate to drop me an E-Mail or open an issue yourself.
I am not responsible for any results produced with VeGETA nor for the conclusions you draw from it.

***
## VeGETA - Viral GEnome sTructure Alignments
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)![Python3.6](https://img.shields.io/badge/Language-Python_3.6-steelblue.svg)[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/klamkiewicz?label=%40klamkiewicz&style=social)](https://twitter.com/klamkiewicz)

***

### What is this about?
RNA viruses usually exploit RNA secondary structures to regulate all kind of mechanisms and steps during the viral life-cycle. For example, [earlier studies](https://www.sciencedirect.com/science/article/pii/S004268221730404X) indicate that a single mutation in the 5'UTR of HCoV-229E, which disrupts a specific structural element, causes the viral replication level to decrease. Therefore we provide this pipeline to generate whole genome alignments of viruses and annotates the different regions with RNA secondary structures information.
### What do I need?
Simply spoken, all you need is a `.fasta` file containing your viruses and a linux OS. As usual in life, things are more complicated. You'll need to install some third party software to use VeGETA. Please refer to the [`How-To Install` section](#how-do-i-install-vegeta).

### What do I get?
VeGETA outputs a whole bunch of alignments. First and foremost, you'll get an alignment with structure information for the representative viruses of your input set. However, most of the time RNA viruses are very diverse in their sequence information. Since we're using a sequence alignment as first scaffold, the diversity usually leads to problems further downstream the pipeline. You will get an alignment - we never promised it'd be beautiful. 

***

### How do I install VeGETA?
The easiest way *at the moment* is simply following these instructions:

Download the latest version of VeGETA or clone this repository:
Download the latest version of VeGETA [here](https://github.com/klamkiew/vegeta/archive/v0.3.zip) or clone this repository:

`git clone https://github.com/klamkiew/vegeta.git`

To make it as easy as possible, I highly recommend using `conda` for the next steps.
[Here](https://docs.conda.io/en/latest/miniconda.html) you can get the conda installer of your choice, depending on the OS you are working with.

Simply create a new environment based on the `environment.yml` file of this repository.
This will create an environment for you called `vegeta` which already includes all dependencies and third-party tools used by VeGETA. Confused? Check our these [conda webpages](https://docs.conda.io/en/latest/) if you want to learn more.

`conda env create -f environment.yml`


Next, activate your newly created environment:

`conda activate vegeta`

And you are all set up: Run a `vegeta.py -h` to see whether any package is missing. Seeing the help message of our pipeline? Good, then we can continue.

### How to use VeGETA
Using VeGETA is relatively easy; simply prepare a `.fasta` file with your sequences and run the following command:

`vegeta <PATH/TO/YOUR/FASTA>`
That's it. Be aware, that this may take some time, depending on how many sequences you want to input. Experience says that you need minutes to some hours for ~100 sequences, several hours for >1.000 sequences. We broke our pipeline with 60.000 sequences, fixes are planned already, but not implemented yet.
To give you some quality of life improvements, we'd advise using the following parameters of VeGETA:

`-p <INT>`: sets the number of used processes. Per default this is 1. 

`-v`: Get some more information about what is going on right now.

`-o <PATH/TO/DIR>`: Specify an output path. Per default, we simply create a new directory called `vegeta/` in your current working directory.

There are some more parameters which are not explained in detail here. Please refer to the help page of VeGETA (`vegeta --help`).

***

### Troubleshooting

**I get a `vegeta: command not found` error**

Well you most likely did not include VeGETA to your `$PATH` variable.
Try `PATH=$PATH:</PATH/TO/VEGETA/>` and eventually put this in your `.bashrc` or `.profile`.

**VeGETA tells me that it can't fold sequences of length 0**

True, because that usually doesn't make much sense. Actually, ViennaRNA is telling you this, when it is called by LocARNA. The problems arises in diverse alignments; since refining the scaffold alignment with structure guidance, we extract fragments of the sequence-based alignment and remove gaps from sub-sequences. Sometimes, one sequence ends up being empty. The "issue" was fixed with the latest LocARNA version (2.0.0RC8): make sure it is installed in your conda environment. This is usually the default version when using conda.

**I receive an `OverflowError` when using VeGETA**

Yup, we are aware of this issue. The problem comes from multiprocessing large amounts of data and is most likely related to Issue [#3](https://github.com/klamkiew/vegeta/issues/3).

***

### Future Features
High Priority
* Fix known bugs
* <s>Implement a nucleotide shuffle for significance [#5](https://github.com/klamkiew/vegeta/issues/5) </s>
* <s>Allow a sequence of interest to be present the whole time [#6](https://github.com/klamkiew/vegeta/issues/6)</s>

Medium Priority
* Improve code documentation
* conda support [#7](https://github.com/klamkiew/vegeta/issues/7)
* docker container [#7](https://github.com/klamkiew/vegeta/issues/7)
* Optimize code and RAM usage

***
