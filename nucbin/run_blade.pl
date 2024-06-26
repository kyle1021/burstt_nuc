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
$fcen1 = 0;
$fcen2 = 0;
$rate = 40;	# sampling rate and bandwidth
$rxgain = 0.;	# dB
@gains = ();
$txgain = 50;	# dB
$name = 'dual';
$channels = '1,2';
$nInp = 2;
$fconf = 'user.conf';
$verbose = 0;
$genh5 = 0;


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
-f1 MHz		set central freq for Rx1
-f2 MHz		set central freq for Rx2
		(both -f1 and -f2 are set to -f or default when omitted)
-r MHz		set the sampling rate and bw in MHz ($rate)
-g GAIN		set the bladeRF gain in dB ($rxgain)
--gains 'G1 G2 ...'   set gains for individual channel
			the number of gains should match the number input channels
			e.g. 4 gains for two devices (each two channels)
--no-meta	use bladeRF-cli instead
--single CH	save only one channel from each device (instead of dual channels)
		CH is either 1 or 2
--conf FILE	specify the user config file ($fconf)

--header	generate the .h5 files after recording finishes

";

die $usage if ($#inp < 0);

while (@inp) {
	$_ = shift @inp;
	if (/^--name$/) {
		$name = shift @inp || die $usage;
	} elsif (/^-f$/) {
		$fcen = shift @inp || die $usage;
	} elsif (/^-f1$/) {
		$fcen1 = shift @inp || die $usage;
	} elsif (/^-f2$/) {
		$fcen2 = shift @inp || die $usage;
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
	} elsif (/^--header$/) {
		$genh5 = 1;
	} else {
		$dur = $_;
	};
};

$now = `date +'%y%m%d_%H%M%S'`; chomp $now;
$fname = "${name}_$now";
$sdir = dirname($fname);    # secondary dir below parent
$sdir =~ s/^.\///;
print "sdir: $sdir\n";

if ($fcen1 == 0) {
	$fcen1 = $fcen;
}
if ($fcen2 == 0) {
	$fcen2 = $fcen;
};


$refin  = "10M";
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
print "debug:: gains: @gains\n";


## regardless of nInp, just setup both channels
## then file config, use only the selected channel
if ($nInp == 2) {
$common = "
set frequency rx1 ${fcen1}M
set frequency rx2 ${fcen2}M
set samplerate rx1 ${rate}M
set samplerate rx2 ${rate}M
set bandwidth rx1 ${rate}M
set bandwidth rx2 ${rate}M

set agc off
";
} elsif ($nInp == 1) {
$common = "
set frequency rx$channels ${fcen}M
set samplerate rx$channels ${rate}M
set bandwidth rx$channels ${rate}M

set agc off
";
};

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
@rdirs = ();
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
		push(@rdirs, $sdir);
	} else {
		$rdir1 = "$parent[$i]/$wd";
                if ($sdir ne '.') {
                    $rdir = "$rdir1/$sdir";
                } else {
                    $rdir = $rdir1;
                };
		push(@rdirs, $rdir1);
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
	if ($nInp == 2) {
		$ch1 = $i*2;
		$ch2 = $ch1 + 1;
		print {$FH[$i]} "set gain rx1 $gains[$ch1]\n";
		print {$FH[$i]} "set gain rx2 $gains[$ch2]\n";
	} elsif ($nInp == 1) {
		$ch1 = $i;
		print {$FH[$i]} "set gain rx$channels $gains[$ch1]\n";
	};
	print {$FH[$i]} "trigger J51-1 rx $rxtrig[$i]\n";
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
	    print {$FH[$i]} "set refin_freq $refin\n";
            print {$FH[$i]} "set clock_ref enable\n";
            print {$FH[$i]} "print refin_freq\n";
        } else {
            print {$FH[$i]} "set clock_ref disable\n";
        };

	print {$FH[$i]} "rx config file='$out[$i]' format=bin n=0 channel=$channels\n";
	print {$FH[$i]} "rx; rx start\n";
};

$now2 = `date +'%y%m%d_%H%M%S'`; chomp $now2;
for $i (0 .. $nDev-1) {
	if ($rxtrig[$i] eq 'master') {
		print {$FH[$i]} "trigger J51-1 rx fire\n";
	};
};
for $i (0 .. $nDev-1) {
	print {$FH[$i]} "echo 'starting: $now2'\n";
};

print "duration: $dur, starting: $now2\n";
sleep($dur);	# wait to collect data

$now3 = `date +'%y%m%d_%H%M%S'`; chomp $now3;
for $i (0 .. $nDev-1) {
	print {$FH[$i]} "echo 'ending: $now3'\n";
	print {$FH[$i]} "rx stop; exit\n";
	close($FH[$i]);
};

if ($genh5) {
	@childs = ();

	$cmd0 = "cli2header.py --verbose";
	for $i (0 .. $nDev-1) {
		$pid = fork();
		if ($pid) { # nonzero $pid for is for the parent
			push (@childs, $pid);
			print "$srv: forked pid: $pid\n";
		} elsif ($pid == 0) {	# the child
			$srv = $server[$i];
			$fmeta = $out[$i];
			$rdir = $rdirs[$i];

			if ($srv eq 'local') {
				$cmd = "$cmd0 $fmeta";
			} else {
				$cmd = "ssh $srv -t \'cd $rdir; ~/miniconda3/bin/python ~/local/bin/$cmd0 $fmeta\'";
			};
			print "$cmd\n";
			`$cmd`;
			exit(0);
		} else {	# something wrong
			die "error forking child: $srv\n";
		};
	};

	foreach (@childs) {
		waitpid($_,0);
		print "$srv: fork finished\n";
	};
};
