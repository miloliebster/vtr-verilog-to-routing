#!/usr/bin/env perl

use Term::ANSIColor;
use File::Basename;

my $VTR_ROOT = $ENV{'VTR_ROOT'};

my @qor_files = (
    'chanx_occupancy.txt',
    'chany_occupancy.txt',
    'packing_pin_util.rpt',
    'report_timing.setup.rpt',
    'report_timing.hold.rpt',
    'report_unconstrained_timing.setup.rpt',
    'report_unconstrained_timing.hold.rpt',
);

@circuit_dirs = glob "$VTR_ROOT/vtr_flow/tasks/dejavu/*/newest/*/*";

`mkdir -p $VTR_ROOT/dejavu_results/diffs`;

my $dirty_diff_count = 0;

foreach my $circuit_dir (@circuit_dirs) {
    foreach my $qor_file (@qor_files) {
        my $circuit = basename($circuit_dir);

        $diff_output = `diff $circuit_dir/common_--memoize_cluster_packings_off/$qor_file $circuit_dir/common_--memoize_cluster_packings_on/$qor_file`;
        open(my $FH, '>', "$VTR_ROOT/dejavu_results/diffs/$circuit.$qor_file.diff") or die "oops\n";
        print $FH $diff_output;
        close($FH);

        next if ($diff_output eq "");

        print "diff of $qor_file for circuit $circuit is dirty!\n";
        $dirty_diff_count++;
    }
}

if ($dirty_diff_count == 0) {
    print color('bold green');
} else {
    print color('bold red');
}
print scalar(@circuit_dirs) - $dirty_diff_count, "/", scalar(@circuit_dirs), " diffs were clean\n";
print color('reset');
