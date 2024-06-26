# burstt_nuc

The bin/ contains tools used in the bladeRF testing environment for the BURSTT project.

The 'kyle_bladeRF_src.tar.gz' is a modified source code of the original bladeRF source (https://github.com/Nuand/bladeRF). Specifically, we use the interactive (and scriptable) command line tool 'bladeRF-cli' to record the raw ADC samples. I have enabled saving the metadata along with the raw ADC samples, which tells us if the samples received are contiguous or if any data have been dropped.
