#!/usr/bin/env perl

use File::Copy;
use List::Util qw(any first);
use Term::ANSIColor;

my $VTR_ROOT = $ENV{'VTR_ROOT'};

sub symlink_newest_runs {
    my @task_dirs = glob "$VTR_ROOT/vtr_flow/tasks/dejavu/*";
    foreach my $task_dir (@task_dirs) {
        my @run_dirs = grep { /$task_dir\/run\d{3}/ } glob "$task_dir/*";
        my @sorted_dirs = map  { $_->[0] }
                          sort { $b->[1] <=> $a->[1] }
                          map  { [ $_, (stat($_))[9] ] }
                          @run_dirs;
        if ($sorted_dirs[0]) {
            `ln -sfT $sorted_dirs[0] $task_dir/newest\n`;
        }
    }
}

sub report_results {
    my %results = ();

    @parse_result_files = glob "$VTR_ROOT/vtr_flow/tasks/dejavu/*/newest/parse_results.txt";
    foreach my $parse_result_file (@parse_result_files) {
        open(my $FH, '<', $parse_result_file);

        for (my $i = 0; my $line = <$FH>; $i++) {
            next if ($i == 0);

            chomp $line;

            @fields = split(/\s/, $line);
            my ($arch)    = $fields[0] =~ /(.*)\.xml/;
            my ($circuit) = $fields[1] =~ /(.*?)(_stratix10_arch_timing)?(\.pre-vpr)?\.blif/;
            my ($memo)    = $fields[2] =~ /common_--memoize_cluster_packings_(.*)/;

            $results{$arch}{$circuit}{$memo} = [ @fields[3..6] ];
        }
    }

    `mkdir -p $VTR_ROOT/dejavu_results`;

    foreach my $arch (keys %results) {
        open(my $FD, '>', "$VTR_ROOT/dejavu_results/$arch.txt");

        my $W0 = 24;
        my $W1 = 12;
        my $W2 = 18;

        my $sum_pack_time_speedup_logs = 0;
        my $sum_vpr_time_speedup_logs  = 0;
        my $sum_pack_mem_increase = 0;
        my $sum_vpr_mem_increase  = 0;

        print $FD "+" . ("-" x $W0) . "+" . ((("-" x ($W1 + $W2 + 1)) . "+") x 4) . "\n";
        print $FD "|"                        . (" " x $W0);
        print $FD "| Packing Runtime (mins)" . (" " x ($W1 + $W2 + 1 - length("Packing Runtime (mins)") - 1));
        print $FD "| VPR Runtime (mins)"     . (" " x ($W1 + $W2 + 1 - length("VPR Runtime (mins)")     - 1));
        print $FD "| Packing Memory (MiB)"   . (" " x ($W1 + $W2 + 1 - length("Packing Memory (MiB)")   - 1));
        print $FD "| VPR Memory (MiB)"       . (" " x ($W1 + $W2 + 1 - length("VPR Memory (MiB)")       - 1));
        print $FD "|\n";

        print $FD "+" . ("-" x $W0) . "+" . ((("-" x $W1) . "+" . ("-" x $W2) . "+") x 4) . "\n";
        print $FD "| Benchmark" . (" " x ($W0 - length("Benchmark") - 1));
        for (my $i = 0; $i < 4; $i++) {
            print $FD "| Baseline" . (" " x ($W1 - length("Baseline") - 1));
            print $FD "| Deja Vu"  . (" " x ($W2 - length("Deja Vu")  - 1));
        }
        print $FD "|\n";

        print $FD "+" . ("=" x $W0) . "+" . ((("=" x $W1) . "+" . ("=" x $W2) . "+") x 4) . "\n";
        foreach my $circuit (sort keys %{$results{$arch}}) {
            my $pack_time_off = $results{$arch}{$circuit}{'off'}->[0] / 60;
            my $pack_time_on  = $results{$arch}{$circuit}{'on' }->[0] / 60;
            my $pack_time_speedup = sprintf("%.1f", $pack_time_off / $pack_time_on);
            $sum_pack_time_speedup_logs += log($pack_time_speedup);
            $pack_time_off = sprintf("%.1f", $pack_time_off);
            $pack_time_on  = sprintf("%.1f", $pack_time_on) . ' (' . (' ' x (4 - length($pack_time_speedup))) . $pack_time_speedup . 'x)';

            my $vpr_time_off  = $results{$arch}{$circuit}{'off'}->[2] / 60;
            my $vpr_time_on   = $results{$arch}{$circuit}{'on' }->[2] / 60;
            my $vpr_time_speedup = sprintf("%.1f", $vpr_time_off / $vpr_time_on);
            $sum_vpr_time_speedup_logs += log($vpr_time_speedup);
            $vpr_time_off = sprintf("%.1f", $vpr_time_off);
            $vpr_time_on  = sprintf("%.1f", $vpr_time_on) . ' (' . (' ' x (4 - length($vpr_time_speedup))) . $vpr_time_speedup . 'x)';

            my $pack_mem_off  = $results{$arch}{$circuit}{'off'}->[1];
            my $pack_mem_on   = $results{$arch}{$circuit}{'on' }->[1];
            my $pack_mem_increase = 100 * (abs($pack_mem_off - $pack_mem_on) / $pack_mem_off);
            $sum_pack_mem_increase += $pack_mem_increase;
            $pack_mem_increase = sprintf("%.0f", $pack_mem_increase);
            $pack_mem_off = sprintf("%.1f", $pack_mem_off);
            $pack_mem_on  = sprintf("%.1f", $pack_mem_on) . ' (+' . (' ' x (3 - length($pack_mem_increase))) . $pack_mem_increase . '%)';

            my $vpr_mem_off  = $results{$arch}{$circuit}{'off'}->[3];
            my $vpr_mem_on   = $results{$arch}{$circuit}{'on' }->[3];
            my $vpr_mem_increase = 100 * (abs($vpr_mem_off - $vpr_mem_on) / $vpr_mem_off);
            $sum_vpr_mem_increase += $vpr_mem_increase; $vpr_mem_increase = sprintf("%.0f", $vpr_mem_increase);
            $vpr_mem_off = sprintf("%.1f", $vpr_mem_off);
            $vpr_mem_on  = sprintf("%.1f", $vpr_mem_on) . ' (+' . (' ' x (3 - length($vpr_mem_increase))) . $vpr_mem_increase . '%)';

            print $FD "| $circuit" . (' ' x ($W0 - length($circuit) - 1));
            print $FD "| " . (' ' x ($W1 - length($pack_time_off) - 2)) . $pack_time_off . ' ';
            print $FD "| " . (' ' x ($W2 - length($pack_time_on)  - 2)) . $pack_time_on  . ' ';
            print $FD "| " . (' ' x ($W1 - length($vpr_time_off)  - 2)) . $vpr_time_off  . ' ';
            print $FD "| " . (' ' x ($W2 - length($vpr_time_on)   - 2)) . $vpr_time_on   . ' ';
            print $FD "| " . (' ' x ($W1 - length($pack_mem_off)  - 2)) . $pack_mem_off  . ' ';
            print $FD "| " . (' ' x ($W2 - length($pack_mem_on)   - 2)) . $pack_mem_on   . ' ';
            print $FD "| " . (' ' x ($W1 - length($vpr_mem_off)   - 2)) . $vpr_mem_off   . ' ';
            print $FD "| " . (' ' x ($W2 - length($vpr_mem_on)    - 2)) . $vpr_mem_on    . ' ';
            print $FD "|\n";
        }

        my $n = scalar keys %{$results{$arch}};
        my $pack_time_speedup_geomean = sprintf("(%4.1fx)", exp($sum_pack_time_speedup_logs / $n));
        my $vpr_time_speedup_geomean  = sprintf("(%4.1fx)", exp($sum_vpr_time_speedup_logs  / $n));
        my $pack_mem_increase_mean    = sprintf("(+%3.0f\%)", $sum_pack_mem_increase / $n);
        my $vpr_mem_increase_mean     = sprintf("(+%3.0f\%)", $sum_vpr_mem_increase  / $n);

        print $FD "+" . ("-" x $W0) . "+" . ((("-" x $W1) . "+" . ("-" x $W2) . "+") x 4) . "\n";
        print $FD "| Geomean" . (' ' x ($W0 - length("Geomean") - 1));
        print $FD "| " . (' ' x ($W1 - 1)) . "| " . (' ' x ($W2 - length($pack_time_speedup_geomean) - 2)) . $pack_time_speedup_geomean . ' ';
        print $FD "| " . (' ' x ($W1 - 1)) . "| " . (' ' x ($W2 - length($vpr_time_speedup_geomean)  - 2)) . $vpr_time_speedup_geomean  . ' ';
        print $FD ("| " . (' ' x ($W1 - 1) . "| " . ' ' x ($W2 - 1))) x 2;
        print $FD "|\n";

        print $FD "+" . ("-" x $W0) . "+" . ((("-" x $W1) . "+" . ("-" x $W2) . "+") x 4) . "\n";
        print $FD "| Mean" . (' ' x ($W0 - length("Mean") - 1));
        print $FD ("| " . (' ' x ($W1 - 1) . "| " . ' ' x ($W2 - 1))) x 2;
        print $FD "| " . (' ' x ($W1 - 1)) . "| " . (' ' x ($W2 - length($pack_mem_increase_mean)    - 2)) . $pack_mem_increase_mean    . ' ';
        print $FD "| " . (' ' x ($W1 - 1)) . "| " . (' ' x ($W2 - length($vpr_mem_increase_mean)     - 2)) . $vpr_mem_increase_mean     . ' ';
        print $FD "|\n";

        print $FD "+" . ("-" x $W0) . "+" . ((("-" x $W1) . "+" . ("-" x $W2) . "+") x 4) . "\n";

        close($FD);
    }
}

sub main {
    my $circuit_size = shift(@ARGV);
    my $num_threads = shift(@ARGV);
    $num_threads = 4 if (not $num_threads);

    my @circuit_size_options = ('XXS', 'XS', 'S', 'M', 'L', 'XL', 'XXL');
    if (not any {$_ eq $circuit_size} @circuit_size_options) {
        print "Usage: launch-dejavu-artifacts.pl CIRCUITSIZE [num_threads]\n\n";
        print "Circuit sizes (ymmv):\n";
        print "    XXS  arm_core on flagship only (for script testing)\n";
        print "    XS   circuits under 20 min\n";
        print "    S    circuits under 100 min\n";
        print "    M    circuits under 200 min\n";
        print "    L    circuits under 1000 min\n";
        print "    XL   circuits under 2000 min\n";
        print "    XXL  circuits up to ~5000 min\n";
        exit;
    }

    $circuit_size_index = 0;
    $circuit_size_index++ while ($circuit_size_options[$circuit_size_index] ne $circuit_size);

    @config_dirs = glob "$VTR_ROOT/vtr_flow/tasks/dejavu/*/config";
    foreach my $config_dir (@config_dirs) {
        open(my $FIN, '<', "$config_dir/config.txt.opt");
        open(my $FOUT, '>', "$config_dir/config.txt");
        while (my $line = <$FIN>) {
            $line =~ s/#XXS// if $circuit_size_index >= 0;
            $line =~ s/#XS//  if $circuit_size_index >= 1;
            $line =~ s/#S//   if $circuit_size_index >= 2;
            $line =~ s/#M//   if $circuit_size_index >= 3;
            $line =~ s/#L//   if $circuit_size_index >= 4;
            $line =~ s/#XL//  if $circuit_size_index >= 5;
            $line =~ s/#XXL// if $circuit_size_index >= 6;
            print $FOUT $line;
        }
        close($FIN);
        close($FOUT);
    }

    my @tasks = (
        'dejavu/verilog_flagship',
    );
    push(@tasks, 'dejavu/koios_flagship') if ($circuit_size_index >= 1);
    push(@tasks, 'dejavu/verilog_7series') if ($circuit_size_index >= 1);
    push(@tasks, 'dejavu/verilog_7series') if ($circuit_size_index >= 1);
    push(@tasks, 'dejavu/titan_titanium') if ($circuit_size_index >= 2);

    my $status = system("$VTR_ROOT/vtr_flow/scripts/run_vtr_task.py @tasks -j $num_threads");
    symlink_newest_runs();

    report_results();

    if (not $status) {
        print "\n";
        print color('bold yellow');
        print "+--------------------------------------------------------------------------------+\n";
        print "| Runtime and memory results can be found in \$VTR_ROOT/dejavu_results            |\n";
        print "|                                                                                |\n";
        print "| To verify diffs of QoR reports are clean, run \$VTR_ROOT/compare-dejavu-qor.pl  |\n";
        print "+--------------------------------------------------------------------------------+\n";
        print color('reset');
    }
}

main();
