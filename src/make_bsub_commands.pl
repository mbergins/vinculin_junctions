#!/usr/bin/perl -w

my @addback_dirs = <../results/Addback_nuclei/*>;
my @commands = &build_commands(@addback_dirs);
foreach (@commands) {
	system($_);
}

my @null_dirs = <../results/LKB_null_nuclei/*>;
@commands = &build_commands(@null_dirs);
foreach (@commands) {
	system($_);
}


sub build_commands {
	my @dirs = @_;
	
	my @commands;
	foreach (@dirs) {
		my $out_file = "$_/run.txt";
		my $error_file = "$_/error.txt";
		my $cmd = "bsub -o $out_file -e $error_file -g /cell_polarity matlab -singleCompThread -nodisplay -r \"process_cell_polarity_data('$_')\"";

		push @commands, $cmd;
	}
	
	return @commands;
}
