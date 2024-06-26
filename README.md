# burstt_nuc

The nucbin/ contains tools used in the bladeRF testing environment for the BURSTT project.

The [kyle_bladeRF_src.tar.gz]([url](https://drive.google.com/file/d/18uhRUc8o6Dia2FSuN7r5pTYUrYXkgiRO/view?usp=sharing)) is a modified source code of the original bladeRF source (https://github.com/Nuand/bladeRF). Specifically, we use the interactive (and scriptable) command line tool 'bladeRF-cli' to record the raw ADC samples. I have enabled saving the metadata along with the raw ADC samples, which tells us if the samples received are contiguous or if any data have been dropped.
