#!/usr/bin/perl
# beehive.pl
# calculate timing of copy gain and clonal expansion
# Input: somatic SNV allele fraction
# Output: Timing estimation with confidence interval of copy gain and clonal expansion as well as estimation of clonality
# Author: Lixing Yang, Ben May Department for Cancer Research, The University of Chicago, Chicago, IL, USA
# Email: lixingyang@uchicago.edu

use strict;
use Getopt::Std;
my $version = '0.1';

my %opts = (b=>0, n=>1000, s=>30);
getopts("i:o:t:b:n:s:p:d:m:h", \%opts);

my $inputfile = $opts{i};
my $output_prefix = $opts{o};
my $eventtype = $opts{t};
my $bootstrap = $opts{b}; # default 0, set to enable bootstrapping
my $repeats = $opts{n}; # number of bootstrap
my $spl_parameter = $opts{s};

# initialize
$output_prefix = $inputfile unless (defined($opts{o}));
my $low_bnd = int(0.05*$repeats)-1;
my $up_bnd = int(0.95*$repeats)-1;
$eventtype =~ /(\d)\+(\d)/;
my ($major_chr, $minor_chr) = ($1, $2);
$spl_parameter = ",df=".$spl_parameter;
my $distance_cutoff = 0.2; # distance of 2 highest peaks, most likely major allele and sub-clonal peaks, need to modify to larger value if sub-clonal peak is small and there is small peak not near major peak
my $distance_cutoff2 = 0.06; # distance of the third and fourth peaks to expected locations
my $distance_cutoff3 = 0.06; # distance of the third or fourth peaks to sub-clonal peak
my $peak_cutoff = 0.2; # minimum intensity of peak
my @pre_peaks = split(",", $opts{p});
@pre_peaks = sort(@pre_peaks);
my @denominators = split(",", $opts{d});
my @multipliers = split(",", $opts{m});
#print "$opts{p}\n@pre_peaks\n";exit;

&print_usage if (defined($opts{h}));
unless (defined($inputfile))
{
	die "Please specify inputfile.\nFor additional help, please run perl beehive.pl -h\n";
}
unless (defined($eventtype) or defined($opts{p}))
{
	die "Please specify event type or define peak locations.\nFor additional help, please run perl beehive.pl -h\n";
}
if (defined($opts{p}))
{
	$peak_cutoff = 0.0000001;
	my $number_pre_peaks = scalar @pre_peaks;
	my $number_denominators = scalar @denominators;
	my $number_multipliers = scalar @multipliers;
	#print "@pre_peaks\t$number_pre_peaks\n@denominators\t$number_denominators\n@multipliers\t$number_multipliers\n";
	if ($number_pre_peaks == 2)
	{
		if ($number_pre_peaks == $number_denominators)
		{
			unless ($number_multipliers == 1)
			{
				die "You have specified $number_pre_peaks peaks, therefore 1 multiplier is required.\nFor additional help, please run perl beehive.pl -h\n";
			}
		}
		else
		{
			die "You have specified $number_pre_peaks peaks, therefore $number_pre_peaks denominators are required.\nFor additional help, please run perl beehive.pl -h\n";
		}
	}
	elsif ($number_pre_peaks == 3)
	{
		if ($number_pre_peaks == $number_denominators)
		{
			unless ($number_multipliers == 3)
			{
				die "You have specified $number_pre_peaks peaks, therefore 3 multipliers are required.\nFor additional help, please run perl beehive.pl -h\n";
			}
		}
		else
		{
			die "You have specified $number_pre_peaks peaks, therefore $number_pre_peaks denominators are required.\nFor additional help, please run perl beehive.pl -h\n";
		}
	}
	elsif ($number_pre_peaks == 4)
	{
		if ($number_pre_peaks == $number_denominators)
		{
			unless ($number_multipliers == 6)
			{
				die "You have specified $number_pre_peaks peaks, therefore 6 multipliers are required.\nFor additional help, please run perl beehive.pl -h\n";
			}
		}
		else
		{
			die "You have specified $number_pre_peaks peaks, therefore $number_pre_peaks denominators are required.\nFor additional help, please run perl beehive.pl -h\n";
		}
	}
	else
	{
		die "Only 2 to 4 pre-defined peaks are allowed.\nFor additional help, please run perl beehive.pl -h\n";
	}
}
die "Input file not found\n" unless (-e $inputfile);

my $time0 = time;
my $newline;
my $bsdir = $output_prefix.'_bootstrap/';
if ($bootstrap)
{
	system "rm -rf $bsdir" if (-e $bsdir);
	system "mkdir $bsdir";
	my @newline;
	open FILE, "<$inputfile";
	while ($newline = <FILE>)
	{
		push @newline, $newline;
	}
	close FILE;
	
	for (my $i=0; $i<$repeats; $i++)
	{
		my $af_bs_file = $bsdir.$i;
		open OUT, ">$af_bs_file";
		foreach (@newline)
		{
			my $pick = int(rand(@newline));
			#print "$pick\n";
			print OUT "$newline[$pick]";
		}
		close OUT;
	}
}

my $rfile = $output_prefix.'.r';
my $pdffile = $output_prefix.'.pdf';
my $splfile = $output_prefix.'.txt';
open OUT, ">$rfile";
print OUT qq~pdf("$pdffile")
par(cex=0.5)
data=read.table("$inputfile")
h=hist(data[,1]/data[,2],breaks=100,col=1,xlim=c(0,1),xaxp=c(0,1,10),xlab=("Allele fraction"),main="")
spl=smooth.spline(h\$mids,h\$counts $spl_parameter)
lines(spl,col=2,lwd=4)
h=hist(data[,1]/data[,2],breaks=100,plot=F)
spl=smooth.spline(h\$mids,h\$density $spl_parameter)
write(spl\$x, append=F, file=\"$splfile\")
write(spl\$y, append=T, file=\"$splfile\")
~;
if ($bootstrap)
{
	for (my $i=0; $i<$repeats; $i++)
	{
		my $bssplfile = $bsdir.$i.'.txt';
		my $af_bs_file = $bsdir.$i;
		my $pdf_bs_file = $bsdir.$i.'.pdf';
		print OUT qq~
		data=read.table("$af_bs_file",sep="\t")
		h=hist(data[,1]/data[,2],breaks=100,plot=F)
		spl=smooth.spline(h\$mids,h\$density $spl_parameter)
		write(spl\$x, append=F, file=\"$bssplfile\")
		write(spl\$y, append=T, file=\"$bssplfile\")
		~;
	}
}
close OUT;
system "Rscript $rfile";
system "rm $rfile";

my ($inbin_count1, $inbin_count2, $inbin_count3, $total_count);
my ($peak_major, $peak_minor, $clonality1, $clonality2, $clonality3, $timing1, $timing2, $timing3);
my (@inbin_count1, @inbin_count2, @inbin_count3, @total_count);
my (@peak_major, @peak_minor, @clonality1, @clonality2, @clonality3, @timing1, @timing2, @timing3);
if ($bootstrap)
{
	if (defined($opts{p}))
	{
		for (my $i=0; $i<$repeats; $i++)
		{
			my $bssplfile = $bsdir.$i.'.txt';
			my $af_bs_file = $bsdir.$i;
			my ($bs_inbin_count1, $bs_inbin_count2, $bs_inbin_count3, $bs_total_count, $bs_timing1, $bs_timing2, $bs_timing3) = &timing($bssplfile, $af_bs_file);
			#print "$i\t$bs_timing1\t$bs_timing2\t$bs_timing3\n";
			push @timing1, $bs_timing1;
			push @timing2, $bs_timing2;
			push @timing3, $bs_timing3;
			@timing1 = sort {$a <=> $b} @timing1;
			@timing2 = sort {$a <=> $b} @timing2;
			@timing3 = sort {$a <=> $b} @timing3;
		}
	}
	else
	{
		for (my $i=0; $i<$repeats; $i++)
		{
			my $bssplfile = $bsdir.$i.'.txt';
			my $af_bs_file = $bsdir.$i;
			my ($bs_inbin_count1, $bs_inbin_count2, $bs_total_count, $bs_peak_major, $bs_peak_minor, $bs_clonality1, $bs_clonality2, $bs_clonality3, $bs_timing1, $bs_timing2, $bs_timing3) = &timing($splfile, $af_bs_file);
			push @timing1, $bs_timing1;
			push @timing2, $bs_timing2;
			push @timing3, $bs_timing3;
			@timing1 = sort {$a <=> $b} @timing1;
			@timing2 = sort {$a <=> $b} @timing2;
			@timing3 = sort {$a <=> $b} @timing3;
		}
	}
}
if (defined($opts{p}))
{
	($inbin_count1, $inbin_count2, $inbin_count3, $total_count, $timing1, $timing2, $timing3) = &timing($splfile, $inputfile);
	print "Number of SNVs in bin1: $inbin_count1\n";
	if ($timing2)
	{
		print "Number of SNVs in bin2: $inbin_count2\n";
	}
	if ($timing3)
	{
		print "Number of SNVs in bin3: $inbin_count3\n";
	}
	print "Total number of SNVs: $total_count\n";
	print "Timing of 1st event: $timing1 (@timing1[$low_bnd], @timing2[$up_bnd])\n";
	if ($timing2)
	{
		print "Timing of 2nd event: $timing2 (@timing2[$low_bnd], @timing1[$up_bnd])\n";
	}
	if ($timing3)
	{
		print "Timing of 3rd event: $timing3 (@timing3[$low_bnd], @timing3[$up_bnd])\n";
	}
}
else
{
	($inbin_count1, $inbin_count2, $total_count, $peak_major, $peak_minor, $clonality1, $clonality2, $clonality3, $timing1, $timing2) = &timing($splfile, $inputfile);
	print "Number of SNVs in bin1: $inbin_count1\n";
	if ($peak_minor)
	{
		print "Number of SNVs in bin2: $inbin_count2\n";
	}
	print "Total number of SNVs: $total_count\n";
	print "Peak position for amplified alleles: $peak_major\n";
	if ($peak_minor)
	{
		print "Peak position for unamplified alleles: $peak_minor\n";
	}
	print "Predicted clonality based on event type: $clonality1\n";
	if ($timing2)
	{
		print "Timing of copy change: $timing1 (@timing1[$low_bnd], @timing2[$up_bnd])\n";
		print "Timing of clonal expansion: $timing2 (@timing2[$low_bnd], @timing2[$up_bnd])\n";
	}
	else
	{
		print "Timing of clonal expansion: $timing1 (@timing1[$low_bnd], @timing1[$up_bnd])\n";
	}
}

system "rm -rf $bsdir" if ($bootstrap);
system "rm $splfile";
my $time1 = time;
my $local1 = localtime($time1);
my $difference = $time1 - $time0;
my $seconds    =  $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes    =  $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "Time used: $difference:$minutes:$seconds\n";


sub timing
{
	my $splfile = shift;
	my $inputfile = shift;
	
	my ($ref_peaks, $ref_cutoffs, $cutoff1, $cutoff2, $cutoff3);
	# cutoff1,2,3 ordered by min->max
	if (defined($opts{p}))
	{
		($ref_peaks, $ref_cutoffs) = &pre_peaks($splfile, \@pre_peaks);
		if (@pre_peaks == 2)
		{
			$cutoff1 = @$ref_cutoffs[0];
			#print "@$ref_cutoffs\ncutoff1: $cutoff1\n";
		}
		if (@pre_peaks == 3)
		{
			($cutoff2, $cutoff1) = @$ref_cutoffs;
			unless ($$ref_peaks[1])
			{
				$cutoff1 = $cutoff2 + 0.01;
				$cutoff2 -= 0.01;
			}
		}
		if (@pre_peaks == 4)
		{
			($cutoff3, $cutoff2, $cutoff1) = @$ref_cutoffs;
			if ($$ref_peaks[1])
			{
				if ($$ref_peaks[2])
				{
					;
				}
				else
				{
					$cutoff1 = $cutoff2 + 0.01;
					$cutoff2 -= 0.01;
				}
			}
			else
			{
				if ($$ref_peaks[2])
				{
					$cutoff1 = $cutoff2;
					$cutoff2 = $cutoff3 + 0.01;
					$cutoff3 -= 0.01;
				}
				else
				{
					$cutoff3 -= 0.02;
					$cutoff2 = $cutoff3 + 0.02;
					$cutoff1 = $cutoff2 + 0.02;
				}
			}
		}
	}
	elsif ($eventtype =~ /1\+0/ or $eventtype =~ /1\+1/)
	{
		# 2 peaks
		($ref_peaks, $ref_cutoffs) = &peaks($splfile, 2);
		$cutoff1 = $$ref_cutoffs[0];
	}
	# X+0, X+1, X+X, X>=2
	elsif ($minor_chr == 0 or $minor_chr == 1 or $major_chr == $minor_chr)
	{
		# 3 peaks
		($ref_peaks, $ref_cutoffs) = &peaks($splfile, 2);
		($ref_peaks, $ref_cutoffs) = &peaks($splfile, 3, $$ref_peaks[1]/$major_chr);
		($cutoff2, $cutoff1) = @$ref_cutoffs;
		unless ($$ref_peaks[1])
		{
			$cutoff1 = $cutoff2 + 0.01;
			$cutoff2 -= 0.01;
		}
	}
	# 3+2, 4+3, 4+2, ...
	else
	{
		# 4 peaks
		($ref_peaks, $ref_cutoffs) = &peaks($splfile, 2);
		($ref_peaks, $ref_cutoffs) = &peaks($splfile, 4, $$ref_peaks[1]/$major_chr, $$ref_peaks[1]*$minor_chr/$major_chr);
		($cutoff3, $cutoff2, $cutoff1) = @$ref_cutoffs;
		if ($$ref_peaks[1])
		{
			if ($$ref_peaks[2])
			{
				;
			}
			else
			{
				$cutoff1 = $cutoff2 + 0.01;
				$cutoff2 -= 0.01;
			}
		}
		else
		{
			if ($$ref_peaks[2])
			{
				$cutoff1 = $cutoff2;
				$cutoff2 = $cutoff3 + 0.01;
				$cutoff3 -= 0.01;
			}
			else
			{
				$cutoff3 -= 0.02;
				$cutoff2 = $cutoff3 + 0.02;
				$cutoff1 = $cutoff2 + 0.02;
			}
		}
	}
	print "Peaks called: @$ref_peaks\nValleys called: $cutoff3 $cutoff2 $cutoff1\n" unless ($inputfile =~ /bootstrap/);
	
	my ($total_count, $inbin_count1, $inbin_count2, $inbin_count3);
	open FILE, "<$inputfile";
	while ($newline = <FILE>)
	{
		chomp $newline;
		my @data = split ("\t", $newline);
		$total_count++;
		my $af = $data[0] / $data[1];
		if ($af >= $cutoff1)
		{
			$inbin_count1++;
		}
		if ($af >= $cutoff2 and $af < $cutoff1)
		{
			$inbin_count2++;
		}
		if ($af >= $cutoff3 and $af < $cutoff2)
		{
			$inbin_count3++;
		}
		
	}
	close FILE;
	
	my ($peak_major, $peak_minor, $peak_subclone, $timing1, $timing2, $timing3, $clonality1, $clonality2, $clonality3);
	# clonality1: major chr, clonality2: minor chr, clonality3: after copy change, before clonal expansion
	
	if (defined($opts{p}))
	{
		my ($inbin_count_adj1, $inbin_count_adj2, $inbin_count_adj3, $inbin_count_adj4, $total_count_adj);
		my $inbin_count4 = $total_count - $inbin_count3 - $inbin_count2 - $inbin_count1;
		$inbin_count_adj1 = $inbin_count1/$denominators[0];
		$inbin_count_adj2 = ($inbin_count2 - $multipliers[0]*$inbin_count_adj1) / $denominators[1];
		$inbin_count_adj3 = ($inbin_count3 - $multipliers[2]*$inbin_count_adj2 - $multipliers[1]*$inbin_count_adj1) / $denominators[2] if (defined($denominators[2]));
		$inbin_count_adj4 = ($inbin_count4 - $multipliers[5]*$inbin_count_adj3 - $multipliers[4]*$inbin_count_adj2 - $multipliers[3]*$inbin_count_adj1) / $denominators[3] if (defined($denominators[3]));
		$total_count_adj = $inbin_count_adj1 + $inbin_count_adj2 + $inbin_count_adj3 + $inbin_count_adj4;
		$timing1 = $inbin_count_adj1 / $total_count_adj; # 1st event
		$timing1 = sprintf("%.2f", $timing1) if (defined($timing1));
		if (defined($pre_peaks[2]))
		{
			$timing2 = ($inbin_count_adj1 + $inbin_count_adj2) / $total_count_adj; # 2nd event
			$timing2 = sprintf("%.2f", $timing2) if (defined($timing2));
		}
		if (defined($pre_peaks[3]))
		{
			$timing3 = ($inbin_count_adj1 + $inbin_count_adj2 + $inbin_count_adj3) / $total_count_adj; # 3rd event
			$timing3 = sprintf("%.2f", $timing3) if (defined($timing3));
		}
		#print "$inbin_count_adj1, $inbin_count_adj2, $inbin_count_adj3, $inbin_count_adj4, $total_count_adj\n";
		return ($inbin_count1, $inbin_count2, $inbin_count3, $total_count, $timing1, $timing2, $timing3);
	}
	else
	{
		if (defined($cutoff1))
		{
			$peak_major = $$ref_peaks[-1];
		}
		if (defined($cutoff2))
		{
			$peak_minor = $$ref_peaks[-2];
		}
		if (defined($cutoff3))
		{
			$peak_subclone = $$ref_peaks[0];
		}
		
		# calculate clonality
		$clonality1 = 2*$peak_major/($major_chr-($major_chr+$minor_chr-2)*$peak_major);
		
		# calculate timing
		if ($minor_chr == 0)
		{
			$timing1 = $inbin_count1 / ($inbin_count1 + ($total_count-$inbin_count1)/$major_chr);
			if (defined($cutoff2))
			{
				my $inbin_count_adj2 = $inbin_count2/$major_chr;
				my $inbin_count_adj3 = ($total_count - $inbin_count2 - $inbin_count1)/$major_chr;
				my $total_count_adj = $inbin_count1 + $inbin_count_adj2 + $inbin_count_adj3;
				$timing1 = $inbin_count1 / $total_count_adj;
				$timing2 = ($inbin_count1 + $inbin_count_adj2) / $total_count_adj;
				$clonality3 = 2*$peak_minor/(1-($major_chr-2)*$peak_minor);
			}
		}
		elsif($minor_chr == 1)
		{
			if ($major_chr == $minor_chr)
			{
				$timing1 = $inbin_count1/2 / ($inbin_count1/2 + ($total_count-$inbin_count1)/($major_chr+$minor_chr));
			}
			else
			{
				$timing2 = $inbin_count1 / ($inbin_count1 + ($total_count-2*$inbin_count1)/($major_chr+$minor_chr));
				if (defined($cutoff2))
				{
					my $inbin_count_adj2 = ($inbin_count2 - $inbin_count1)/($major_chr+$minor_chr);
					my $inbin_count_adj3 = ($total_count - $inbin_count2 - $inbin_count1)/($major_chr+$minor_chr);
					my $total_count_adj = $inbin_count1 + $inbin_count_adj2 + $inbin_count_adj3;
					$timing1 = $inbin_count1 / $total_count_adj;
					$timing2 = ($inbin_count1 + $inbin_count_adj2) / $total_count_adj;
					$clonality2 = 2*$peak_minor/($minor_chr-($major_chr+$minor_chr-2)*$peak_minor);
				}
			}
		}
		else
		{
			if ($major_chr == $minor_chr)
			{
				$timing1 = $inbin_count1/2 / ($inbin_count1/2 + ($total_count-$inbin_count1)/($major_chr+$minor_chr));
				if (defined($cutoff2))
				{
					my $inbin_count_adj2 = $inbin_count2/($major_chr+$minor_chr);
					my $inbin_count_adj3 = ($total_count - $inbin_count2 - $inbin_count1)/($major_chr+$minor_chr);
					my $total_count_adj = $inbin_count1/2 + $inbin_count_adj2 + $inbin_count_adj3;
					$timing1 = $inbin_count1/2 / $total_count_adj;
					$timing2 = ($inbin_count1/2 + $inbin_count_adj2) / $total_count_adj;
				}
			}
			else
			{
				$timing1 = $inbin_count1 / ($inbin_count1 + ($total_count-2*$inbin_count1)/($major_chr+$minor_chr));
				if (defined($cutoff2) and defined($cutoff3))
				{
					my $inbin_count_adj3 = $inbin_count3 / ($major_chr+$minor_chr);
					my $inbin_count_adj4 = ($total_count - $inbin_count3 - $inbin_count2 - $inbin_count1) / ($major_chr+$minor_chr);
					my $total_count_adj = ($inbin_count1 + $inbin_count2)/2 + $inbin_count_adj3 + $inbin_count_adj4;
					$timing1 = ($inbin_count1 + $inbin_count2)/2 / $total_count_adj;
					$timing2 = (($inbin_count1 + $inbin_count2)/2 + $inbin_count_adj3) / $total_count_adj;
					$clonality2 = 2*$peak_minor/($minor_chr-($major_chr+$minor_chr-2)*$peak_minor);
					$clonality3 = 2*$peak_subclone/(1-($major_chr+$minor_chr-2)*$peak_subclone);
				}
			}
		}
		
		$timing1 = sprintf("%.2f", $timing1) if (defined($timing1));
		$timing2 = sprintf("%.2f", $timing2) if (defined($timing2));
		$clonality1 = sprintf("%.2f", $clonality1) if (defined($clonality1));
		$clonality2 = sprintf("%.2f", $clonality2) if (defined($clonality2));
		$clonality3 = sprintf("%.2f", $clonality3) if (defined($clonality3));
		
		return ($inbin_count1, $inbin_count2, $total_count, $peak_major, $peak_minor, $clonality1, $clonality2, $clonality3, $timing1, $timing2);
	}
}

sub pre_peaks
# call peaks on pre-defined locations
{
	my $splfile = shift;
	my $ref_pre_peaks = shift;
	
	my (@alldata, @mid, @intensity, %all_peaks, %all_valleys, $highest_peak, $highest_intensity, $second_high_peak, $second_high_intensity);
	open FILE, "<$splfile";
	while ($newline = <FILE>)
	{
		chomp $newline;
		my @data = split (" ", $newline);
		push @alldata, @data;
	}
	close FILE;
	for (my $i=0;$i<@alldata/2;$i++)
	{
		push @mid, $alldata[$i];
	}
	for (my $i=@alldata/2;$i<@alldata;$i++)
	{
		push @intensity, $alldata[$i];
	}
	if (length($mid[0]) == 5)
	{
		foreach (@mid)
		{
			$_ += 0.005; # round up
		}
	}
	
	for (my $i=1; $intensity[$i]; $i++)
	{
		next if ($intensity[$i] < $peak_cutoff);
		if ($intensity[$i] > $intensity[$i-1] and $intensity[$i] > $intensity[$i+1])
		{
			#print "peak: $mid[$i]\t$intensity[$i]\n";
			$all_peaks{$mid[$i]}[0] = $intensity[$i];
		}
	}
	for (my $i=1; $intensity[$i]; $i++)
	{
		if ($intensity[$i] < $intensity[$i-1] and $intensity[$i] < $intensity[$i+1])
		{
			#print "valley: $mid[$i]\t$intensity[$i]\n";
			$all_valleys{$mid[$i]}[0] = $intensity[$i];
		}
	}
	
	# %all_peaks, %all_valleys, key: position of peak
	# 0, intensity of peak
	my (@peaks, @cutoffs);
	foreach (@$ref_pre_peaks)
	{
		my $pre_peak = $_;
		my $called_peak;
		foreach my $peak (keys %all_peaks)
		{
			if (abs($peak - $pre_peak) < abs($called_peak - $pre_peak))
			{
				$called_peak = $peak;
			}
		}
		if (abs($called_peak - $pre_peak) > $distance_cutoff2)
		{
			$called_peak = '';
		}
		push @peaks, $called_peak;
		#print "defined peak: $pre_peak\tcalled peak: $called_peak\n";
	}
		
	# find valleys between called peaks
	for (my $i=0;($peaks[$i+1] or $peaks[$i+2] or $peaks[$i+3]);$i++)
	{
		my $lowest_valley;
		my $lowest_valley_intensity = 100;
		#print "$i\t$peaks[$i]\t$peaks[$i+1]\t$peaks[$i+2]\tyes\n";
		next unless ($peaks[$i]);
		if ($peaks[$i+1])
		{
			foreach my $valley (keys %all_valleys)
			{
				if ($valley > $peaks[$i] and $valley < $peaks[$i+1] and $all_valleys{$valley}[0] < $lowest_valley_intensity)
				{
					#print "$valley\t$lowest_valley\tyes\n";
					$lowest_valley_intensity = $all_valleys{$valley}[0];
					$lowest_valley = $valley;
				}
			}
		}
		elsif ($peaks[$i+2])
		{
			foreach my $valley (keys %all_valleys)
			{
				if ($valley > $peaks[$i] and $valley < $peaks[$i+2] and $all_valleys{$valley}[0] < $lowest_valley_intensity)
				{
					$lowest_valley_intensity = $all_valleys{$valley}[0];
					$lowest_valley = $valley;
				}
			}
		}
		elsif ($peaks[$i+3])
		{
			foreach my $valley (keys %all_valleys)
			{
				if ($valley > $peaks[$i] and $valley < $peaks[$i+3] and $all_valleys{$valley}[0] < $lowest_valley_intensity)
				{
					$lowest_valley_intensity = $all_valleys{$valley}[0];
					$lowest_valley = $valley;
				}
			}
		}
		push @cutoffs, $lowest_valley;
	}
	#print "@peaks\n@cutoffs\n";
	return (\@peaks, \@cutoffs);
}

sub peaks
# call 2 major peaks first, then call 3rd and 4th peak based on pre-defined locations
# if 2 peaks, 0 sub-clonal peak, 1 clonal peak
# if 3 peaks, 0 sub-clonal peak, 1 minor allele peak, 2 major allele peak
# if 4 peaks, 0 sub-clonal peak, 1 single allele peak, 2 minor allele peak, 3 major allele peak
{
	my $splfile = shift;
	my $number_of_peaks = shift;
	my $third_peak_location = shift;
	my $fourth_peak_location = shift;
	#print "$number_of_peaks\n";
	
	my (@alldata, @mid, @intensity, %all_peaks, %all_valleys, $highest_peak, $highest_intensity, $second_high_peak, $second_high_intensity);
	open FILE, "<$splfile";
	while ($newline = <FILE>)
	{
		chomp $newline;
		my @data = split (" ", $newline);
		push @alldata, @data;
	}
	close FILE;
	for (my $i=0;$i<@alldata/2;$i++)
	{
		push @mid, $alldata[$i];
	}
	for (my $i=@alldata/2;$i<@alldata;$i++)
	{
		push @intensity, $alldata[$i];
	}
	# round up
	foreach (@mid)
	{
		if (length($_) == 5)
		{
				$_ += 0.005;
		}
		else
		{
			$_ = sprintf("%.2f", $_);
		}
	}
	
	$highest_peak = $mid[0];
	$highest_intensity = $intensity[0];
	for (my $i=1; $intensity[$i]; $i++)
	{
		next if ($intensity[$i] < $peak_cutoff);
		if ($intensity[$i] > $intensity[$i-1] and $intensity[$i] > $intensity[$i+1])
		{
			#print "peak: $mid[$i]\t$intensity[$i]\n";
			$all_peaks{$mid[$i]}[0] = $intensity[$i];
		}
		if ($intensity[$i] > $highest_intensity)
		{
			$highest_peak = $mid[$i];
			$highest_intensity = $intensity[$i];
		}
	}
	for (my $i=1; $intensity[$i]; $i++)
	{
		if ($intensity[$i] < $intensity[$i-1] and $intensity[$i] < $intensity[$i+1])
		{
			#print "valley: $mid[$i]\t$intensity[$i]\n";
			$all_valleys{$mid[$i]}[0] = $intensity[$i];
		}
	}
	
	# %all_peaks, %all_valleys, key: position of peak
	# 0, intensity of peak
	# 1, position of left peak
	# 2, position of right peak
	# 3, distance to highest peak
	# 4, distance to pre-defined third peak location
	# 5, distance to pre-defined fourth peak location
	foreach my $peak (keys %all_peaks)
	{
		my ($left, $right, $left_distance, $right_distance);
		foreach my $valley (keys %all_valleys)
		{
			if ($valley < $peak)
			{
				if ($left_distance)
				{
					if ($peak - $valley < $left_distance)
					{
						$left_distance = $peak - $valley;
						$left = $valley;
					}
				}
				else
				{
					$left_distance = $peak - $valley;
					$left = $valley;
				}
				$all_peaks{$peak}[1] = $left;
			}
			if ($valley > $peak)
			{
				if ($right_distance)
				{
					if ($valley - $peak < $right_distance)
					{
						$right_distance = $valley - $peak;
						$right = $valley;
					}
				}
				else
				{
					$right_distance = $valley - $peak;
					$right = $valley;
				}
				$all_peaks{$peak}[2] = $right;
			}
		}
		$all_peaks{$peak}[3] = abs($highest_peak - $peak);
		$all_peaks{$peak}[4] = abs($third_peak_location - $peak);
		$all_peaks{$peak}[5] = abs($fourth_peak_location - $peak);
		if ($peak ne $highest_peak)
		{
			if ($all_peaks{$peak}[3] >= $distance_cutoff)
			{
				if ($all_peaks{$peak}[0] > $second_high_intensity)
				{
					$second_high_intensity = $all_peaks{$peak}[0];
					$second_high_peak = $peak;
				}
			}
		}
		#print "peak: $peak\t$all_peaks{$peak}[1]\t$all_peaks{$peak}[2]\t$all_peaks{$peak}[3]\n";
	}
	foreach my $valley (keys %all_valleys)
	{
		my ($left, $right, $left_distance, $right_distance);
		foreach my $peak (keys %all_peaks)
		{
			if ($peak < $valley)
			{
				if ($left_distance)
				{
					if ($valley - $peak < $left_distance)
					{
						$left_distance = $valley - $peak;
						$left = $peak;
					}
				}
				else
				{
					$left_distance = $valley - $peak;
					$left = $peak;
				}
				$all_valleys{$valley}[1] = $left;
			}
			if ($peak > $valley)
			{
				if ($right_distance)
				{
					if ($peak - $valley < $right_distance)
					{
						$right_distance = $peak - $valley;
						$right = $peak;
					}
				}
				else
				{
					$right_distance = $peak - $valley;
					$right = $peak;
				}
				$all_valleys{$valley}[2] = $right;
			}
		}
		#print "valley: $valley\t$all_valleys{$valley}[1]\t$all_valleys{$valley}[2]\n";
	}
	#print "Highest peak: $highest_peak\nSecond highest peak: $second_high_peak\n";
	
	my (@peaks, @cutoffs);
	if ($number_of_peaks >= 2)
	{
		@peaks = ($highest_peak, $second_high_peak);
		@peaks = sort {$a <=> $b} @peaks;
	}
	my ($third_peak, $fourth_peak, $third_peak_intensity, $fourth_peak_intensity);
	my ($third_peak_distance, $fourth_peak_distance) = (1, 1);
	if ($number_of_peaks >= 3 and defined($third_peak_location))
	{
		# the peak with given location is already called
		if (abs($highest_peak - $third_peak_location) <= $distance_cutoff2 or abs($second_high_peak - $third_peak_location) <= $distance_cutoff2)
		{
			foreach my $peak (keys %all_peaks)
			{
				if ($peak ne $highest_peak and $peak ne $second_high_peak)
				{
					if (abs($peak - $highest_peak) > $distance_cutoff2 and abs($peak - $second_high_peak) > $distance_cutoff2)
					{
						if ($all_peaks{$peak}[0] > $third_peak_intensity and $peak < $highest_peak and $peak < $second_high_peak)
						{
							$third_peak_intensity = $all_peaks{$peak}[0];
							$third_peak = $peak;
						}
					}
				}
			}
		}
		else
		{
			foreach my $peak (keys %all_peaks)
			{
				if ($peak ne $highest_peak and $peak ne $second_high_peak)
				{
					if ($all_peaks{$peak}[4] <= $third_peak_distance and abs($peak - $third_peak_location) <= $distance_cutoff2)
					{
						$third_peak_distance = $all_peaks{$peak}[4];
						$third_peak = $peak;
					}
				}
			}
		}
		# minor allele peak overlaps sub-clonal peak
		unless ($third_peak)
		{
				if (abs($peaks[0] - $third_peak_location) <= $distance_cutoff3)
			{
				$third_peak = $peaks[0];
			}
		}
		#print "third peak: $third_peak\n";
		push @peaks, $third_peak;
		@peaks = sort {$a <=> $b} @peaks;
		
		# peak too small to call
		unless ($peaks[0])
		{
			($peaks[0], $peaks[1]) = ($peaks[1], $peaks[0]);
		}
		#print "$peaks[0]\t$peaks[1]\t$peaks[2]\n";
	}
	if ($number_of_peaks >= 4 and defined($third_peak_location) and defined($fourth_peak_location))
	{
		# the peak with given location is already called
		if (abs($highest_peak - $fourth_peak_location) <= $distance_cutoff2 or abs($second_high_peak - $fourth_peak_location) <= $distance_cutoff2 or abs($third_peak - $fourth_peak_location) <= $distance_cutoff2)
		{
			foreach my $peak (keys %all_peaks)
			{
				if ($peak ne $highest_peak and $peak ne $second_high_peak and $peak ne $third_peak)
				{
					if (abs($peak - $highest_peak) >= $distance_cutoff and abs($peak - $second_high_peak) >= $distance_cutoff and abs($peak - $third_peak) >= $distance_cutoff)
					{
						if ($all_peaks{$peak}[0] > $fourth_peak_intensity and $peak < $highest_peak and $peak < $second_high_peak and $peak > $third_peak)
						{
							$fourth_peak_intensity = $all_peaks{$peak}[0];
							$fourth_peak = $peak;
						}
					}
				}
			}
		}
		else
		{
			foreach my $peak (keys %all_peaks)
			{
				if ($peak ne $highest_peak and $peak ne $second_high_peak and $peak ne $third_peak)
				{
					if ($all_peaks{$peak}[5] < $fourth_peak_distance)
					{
						$fourth_peak_distance = $all_peaks{$peak}[5];
						$fourth_peak = $peak;
					}
				}
			}
		}
		push @peaks, $fourth_peak;
		@peaks = sort {$a <=> $b} @peaks;
		
		# peak too small to call
		unless ($peaks[0])
		{
			if ($peaks[1])
			{
				if (abs($peaks[2] - $third_peak_location) < abs($peaks[2] - $fourth_peak_location))
				{
					($peaks[0], $peaks[1], $peaks[2]) = ($peaks[1], $peaks[2], $peaks[0]);
				}
				else
				{
					($peaks[0], $peaks[1]) = ($peaks[1], $peaks[0]);
				}
			}
			# third and fourth peaks both too small to call
			else
			{
				($peaks[0], $peaks[2]) = ($peaks[2], $peaks[0]);
			}
		}
	}
	
	# find valleys between called peaks
	for (my $i=0;($peaks[$i+1] or $peaks[$i+2] or $peaks[$i+3]);$i++)
	{
		my $lowest_valley;
		my $lowest_valley_intensity = 100;
		#print "$i\t$peaks[$i]\t$peaks[$i+1]\t$peaks[$i+2]\tyes\n";
		next unless ($peaks[$i]);
		if ($peaks[$i+1])
		{
			foreach my $valley (keys %all_valleys)
			{
				if ($valley > $peaks[$i] and $valley < $peaks[$i+1] and $all_valleys{$valley}[0] < $lowest_valley_intensity)
				{
					$lowest_valley_intensity = $all_valleys{$valley}[0];
					$lowest_valley = $valley;
				}
			}
		}
		elsif ($peaks[$i+2])
		{
			foreach my $valley (keys %all_valleys)
			{
				if ($valley > $peaks[$i] and $valley < $peaks[$i+2] and $all_valleys{$valley}[0] < $lowest_valley_intensity)
				{
					$lowest_valley_intensity = $all_valleys{$valley}[0];
					$lowest_valley = $valley;
				}
			}
		}
		elsif ($peaks[$i+3])
		{
			foreach my $valley (keys %all_valleys)
			{
				if ($valley > $peaks[$i] and $valley < $peaks[$i+3] and $all_valleys{$valley}[0] < $lowest_valley_intensity)
				{
					$lowest_valley_intensity = $all_valleys{$valley}[0];
					$lowest_valley = $valley;
				}
			}
		}
		push @cutoffs, $lowest_valley;
	}
	if ($number_of_peaks >= 3)
	{
		# minor allele peak overlaps with sub-clonal peak, set cutoff to left most position
		if ($peaks[0] == $peaks[1])
		{
			$cutoffs[0] = $mid[0];
			#print "cutoff $cutoffs[0], $cutoffs[1]\n";
		}
	}
	#print "@peaks\n@cutoffs\n";
	return (\@peaks, \@cutoffs);
}

sub print_usage
{
	die "beehive.pl [options]
	-i FILE	input file, required
	-o FILE	output file prefix, default input file name
	-t STR	event type, format x+y, x: copy number of major allele, y: copy number for minor allele
	-b INT	bootstrap, 0/1, default 0, set to 1 to enable bootstrap
	-n INT	number of bootstrap, default 1000
	-s INT	spline smooth parameter, default 30
	-p STR	predefined peak locations, minimum 2 peaks, maximum 4 peaks, format: p1,p2,p3,p4
	-d STR	denominators for calculating timing based on predefined peaks, minimum 2, maximum 4, format: d1,d2,d3,d4
	-m STR	multipliers for calculating timing based on predefined peaks, minimum 1, maximum 6, format: m11,m21,m22,m31,m32,m33
	-h help
	version: $version
";
}
