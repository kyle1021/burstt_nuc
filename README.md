# burstt_nuc

The nucbin/ contains tools used in the bladeRF testing environment for the BURSTT project.

The [kyle_bladeRF_src.tar.gz](https://drive.google.com/file/d/18uhRUc8o6Dia2FSuN7r5pTYUrYXkgiRO/view?usp=sharing) is a modified source code of the original bladeRF source (https://github.com/Nuand/bladeRF). Specifically, we use the interactive (and scriptable) command line tool 'bladeRF-cli' to record the raw ADC samples. I have enabled saving the metadata along with the raw ADC samples, which tells us if the samples received are contiguous or if any data have been dropped.

## Preparation of NUC
A [list of components](https://docs.google.com/spreadsheets/d/1VdITM-I1rzFwJ6aPvvJIadkc7IJpfaj7l9BfUYO-Ocw/edit?usp=drive_link) used in the prototype SDR system.
[Connection Schematic](https://docs.google.com/presentation/d/1b8Fiw5of81N5txna7kH7gqhRHJBrmltjTdzQD9ONgUU/edit?usp=drive_link) and some [notes on the set up of the NUC computer](https://docs.google.com/document/d/1KqSft71pUV-_kRTfO7QsiWVZiwI3MP6r7PdH7ggc2nM/edit?usp=drive_link)

## Installation
1. The codes were used in an Ubuntu Linux system (20.04). It has not been tested in any other OS.
2. Download and extract the 'kyle_bladeRF_src.tar.gz' into the home directory. Compile the codes in bladeRF/host/build as instructed in the https://github.com/Nuand/bladeRF
3. clone the nucbin as $HOME/local/bin

## Using the prototype to measure Sun fringes
[Notes on how to take data and use the analysis script](https://docs.google.com/document/d/1ZIqCqajnW9AG7PWrJPK1O2QHiU7ce1K6tyQ9E6kiiZE/edit?usp=drive_link).
