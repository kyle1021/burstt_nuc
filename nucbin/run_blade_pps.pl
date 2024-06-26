#!/usr/bin/perl -w
use Env;
use File::Basename;
use warnings;
no warnings 'once';

# sends master clock from one bladeRF
# the other bladeRF uses it as the external clock
# delay the 2nd bladeRF by 3s before using the clock


#-- defaults --
$clibin = "~/local/bin/bladeRF-meta";	# either 'bladeRF-meta' or 'bladeRF-cli'
$ext  = 'meta';
$fcen = 560;	# TV frequency
$rate = 40;	# sampling rate and bandwidth
$rxgain = 20.;	# dB
@gains = ();
$txgain = 50;	# dB
$name = 'dual';
$channels = '1,2';
$nInp = 2;
$fconf = 'user.conf';
$verbose = 0;


#-- usage --
@inp = @ARGV;
$pg = $0;
$usage = "
record data with multiple bladeRF devices across different servers
enable clock_out on the master board
delay 3s and then 
set clock_sel as external for the slave board

the server and device information should be provided in a file: 'user.conf'

syntax:
$pg <DUR> [options]

options are:
--name NAME	replace the default name, omit the .bin ($name)
-f MHz		set the central freq in MHz ($fcen)
-r MHz		set the sampling rate and bw in MHz ($rate)
-g GAIN		set the bladeRF gain in dB ($rxgain)
--gains 'G1 G2 ...'   set gains for individual channel
			the number of gains should match the number input channels
			e.g. 4 gains for two devices (each two channels)
--no-meta	use bladeRF-cli instead
--single CH	save only one channel from each device (instead of dual channels)
		CH is either 1 or 2
--conf FILE	specify the user config file ($fconf)

";

die $usage if ($#inp < 0);

while (@inp) {
	$_ = shift @inp;
	if (/^--name$/) {
		$name = shift @inp || die $usage;
	} elsif (/^-f$/) {
		$fcen = shift @inp || die $usage;
	} elsif (/^-r$/) {
		$rate = shift @inp || die $usage;
	} elsif (/^-g$/) {
		$rxgain = shift @inp || die $usage;
        } elsif (/^--gains$/) {
                @gains = split /\s+/, shift @inp || die $usage;
	} elsif (/^--no-meta$/) {
		$clibin = '/usr/bin/bladeRF-cli';
		$ext = 'bin';
	} elsif (/^--single$/) {
		$channels = shift @inp || die $usage;
		$nInp = 1;
		if ($name eq 'dual') {
			$name = 'single';
		};
	} elsif (/^--conf$/) {
		$fconf = shift @inp || die $usage;
	} else {
		$dur = $_;
	};
};

$now = `date +'%y%m%d_%H%M%S'`; chomp $now;
$fname = "${name}_$now";
$sdir = dirname($fname);    # secondary dir below parent
$sdir =~ s/^.\///;
print "sdir: $sdir\n";


@clkref = (0);          # backward compatible with older user.conf
if (-s "./$fconf") {
	require "./$fconf";	# the 'user.conf' needs to tbe where the script is executed
} else {
	die "config file not found: $fconf\n";
};


if ($verbose) {
	print "known devices:\n";
   	for (keys %serial) {
        	print "$_ --> $serial{$_}\n";
    	};
	print "\n";
};
$nDev = $#active + 1;	# num of active devices
if ($verbose) {
	print "connections:\n";
	for $i (0 .. $nDev-1) {
		$dev = $active[$i];
		$srv = $server[$i];
		print "$dev --> $srv\n";
	};
	print "\n";
};
if ($#clkref < $#active) {  # backward compatible with older user.conf
    for $i ($#clkref .. $#active-1) {
        push @clkref, 0;
    };
};


$nTotInp = $nDev*$nInp;
if ($#gains < $nTotInp-1) {  # fill @gains with $rxgain if not enough given
	for $i ($#gains+1 .. $nTotInp-1) {
		push @gains, $rxgain;
	};
};


## regardless of nInp, just setup both channels
## then file config, use only the selected channel
$common = "
set frequency rx1 ${fcen}M
set frequency rx2 ${fcen}M
set samplerate rx1 ${rate}M
set samplerate rx2 ${rate}M
set bandwidth rx1 ${rate}M
set bandwidth rx2 ${rate}M

set agc off
";
## remove the gains from common
#set gain rx1 $rxgain
#set gain rx2 $rxgain

## TX worsk with bladeRF-cli (not bladeRF-meta)
$txcen = $fcen + 0.1;
$txonly = "
set frequency tx1 ${txcen}M
set gain tx1 $txgain
tx config file='master_tx.csv' format=csv channel=1 repeat=10
";


@log = ();
@out = ();
@FH  = ();
for $i (0 .. $nDev-1) {
	$dev = $active[$i];
	$ser = $serial{$dev};
	$ifname = "${fname}_$dev";
	$ilog = "$ifname.log";
	push @log, "$ilog";
	push @out, "$ifname.$ext";

	$srv = $server[$i];
	$cmd0 = "$clibin -d \"*:Serial=$ser\"";
	if ($srv eq 'local') {
		$cmd = "$cmd0";	# local, and already in $wd
		$rquote = '';
	} else {
		$rdir1 = "$parent[$i]/$wd";
                if ($sdir ne '.') {
                    $rdir = "$rdir1/$sdir";
                } else {
                    $rdir = $rdir1;
                };
                #print "remote: rdir $rdir\n";
		$cmd = "mkdir -p $rdir; cd $rdir1;";
		$cmd = "$cmd $cmd0";
		$cmd = "ssh $srv -t \'$cmd";
		$rquote = "\'";
	};

	#$FH[$i] = "|$cmd -i >> $ilog\'";	# testing
	open $FH[$i], "|$cmd -i > $ilog 2>& 1$rquote";
	$FH[$i]->autoflush(1);

	print {$FH[$i]} "echo 'Serial: $ser'\n";
	print {$FH[$i]} "echo 'Duration: $dur sec'\n";
	print {$FH[$i]} $common;
	$ch1 = $i*2;
	$ch2 = $ch1 + 1;
	print {$FH[$i]} "set gain rx1 $gains[$ch1]\n";
	print {$FH[$i]} "set gain rx2 $gains[$ch2]\n";
	$cout = $clkout[$i];
	if ($cout) {	# == 1
		print {$FH[$i]} "set clock_out on\n";
	};
};
#print "@log, @out, @FH\n";

sleep(3);	# wait to ensure clock has been sent out

for $i (0 .. $nDev-1) {
	$csel = $clksel[$i];
        print {$FH[$i]} "set clock_sel $csel\n";

        $cref = $clkref[$i]; 
        if ($cref == 1) {
            print {$FH[$i]} "set clock_ref enable\n";
            print {$FH[$i]} "print refin_freq\n";
        } else {
            print {$FH[$i]} "set clock_ref disable\n";
        };

	print {$FH[$i]} "rx config file='$out[$i]' format=bin n=0 channel=$channels timeout=3000\n";
	print {$FH[$i]} "rx\n";
};

# as short as possible commands to arm the triggers
for $i (0 .. $nDev-1) {
	print {$FH[$i]} "trigger J51-1 rx $rxtrig[$i]\n";
	print {$FH[$i]} "rx start\n";
};
## waiting for PPS or trigger

sleep(1);

$now2 = `date +'%y%m%d_%H%M%S'`; chomp $now2;
for $i (0 .. $nDev-1) {
	if ($rxtrig[$i] eq 'master') {
		print {$FH[$i]} "trigger J51-1 rx fire\n";
	};
};
for $i (0 .. $nDev-1) {
	print {$FH[$i]} "echo 'starting: $now2'\n";
};

sleep($dur);	# wait to collect data

$now3 = `date +'%y%m%d_%H%M%S'`; chomp $now3;
for $i (0 .. $nDev-1) {
	print {$FH[$i]} "echo 'ending: $now3'\n";
	print {$FH[$i]} "rx stop; exit\n";
	close($FH[$i]);
};

