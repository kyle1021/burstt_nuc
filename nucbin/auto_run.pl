#!/usr/bin/perl

#$name = "48VD2A";
$name = "LabAC";
$f0 = 190;  # starting central freq
$sep = 20;  # separation between central freq
$g0 = 10;   # gain in dB

for $i (0 .. 8) {
    $f = $f0 + $i*$sep;
    print "$i, $f\n";


    # take 20s data with rate=20Msps
    `./run_blade_meta-nuc-a4.pl 20 -r 20 -f $f -g $g0 --name "${name}_${f}MHz"`;

};

`quick_look.py $name*.meta`;

`plot_waterfall2.py $name*.meta`;


