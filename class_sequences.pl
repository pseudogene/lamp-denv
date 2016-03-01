#!/usr/bin/perl
# $Revision: 0.5 $
# $Date: 2016/03/01 $
# $Id: class_sequences.pl $
# $Author: Michael Bekaert $
#
#
# Copyright 2014-2016, Michaël Bekaert <michael.bekaert@stir.ac.uk>
#
# This file is part of lamp-denv.
#
# lamp-denv is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lamp-denv is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with lamp-denvs. If not, see <http://www.gnu.org/licenses/>.
#
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw/uniq/;
use Bio::SeqIO;
use File::Temp qw/ tempfile /;

#----------------------------------------------------------
our ($VERSION) = 0.5;
my $pathlava = '/lava-dna';

#----------------------------------------------------------
my ($verbose, $combi, $loose, $fasta, $list) = (0, 0, 0);
GetOptions('a|align|alignment=s' => \$fasta, 'l|list=s' => \$list, 'c!' => \$combi, 'e|extra!' => \$loose, 'v|verbose!' => \$verbose);
if (defined $list && -r $list && defined $fasta && -r $fasta)
{
    my (%pop, %sample);
    if (open my $IN, '<', $list)
    {
        while (<$IN>)
        {
            chomp;
            my @data = split m/\t/;
            if (scalar @data >= 6 && defined $data[0] && length($data[0]) > 5 && defined $data[4] && int($data[4]) > 0 && defined $data[5] && int($data[5]) > 0)
            {
                $data[0] = $1 if ($data[0] =~ m/(^\w+\d+\|\d{4})/);
                $sample{$data[0]} = $data[($loose ? 4 : 5)];
                if   (exists $pop{$data[($loose ? 4 : 5)]}) { $pop{$data[($loose ? 4 : 5)]}++; }
                else                                        { $pop{$data[($loose ? 4 : 5)]} = 1; }
            }
        }
        close $IN;
    }
    my $seq_in = Bio::SeqIO->new(-file => '<' . $fasta, -format => 'fasta');
    my %seq_out;
    foreach my $i (keys %pop) { $seq_out{$i} = Bio::SeqIO->new(-file => '>' . $fasta . q{.} . $i . '.fasta', -format => 'fasta'); }
    while (my $inseq = $seq_in->next_seq)
    {
        my $id = $inseq->id();
        $id = $1 if ($id =~ m/(^\w+\d+\|\d{4})/);
        if   (exists $sample{$id}) { $seq_out{$sample{$id}}->write_seq($inseq); }
        else                       { print {*STDERR} '! Unknown sequence ', $inseq->id(), "\n"; }
    }
    foreach my $i (keys %pop) { $seq_out{$i}->close(); }
    $seq_in->close();
    if (scalar keys %pop > 0 && $combi)
    {
        #test combinaisons..
        my @bad;
        my @list = sort keys %pop;
        my $size = scalar @list;
        for (my $i = 0 ; $i < 2**$size ; $i++)
        {
            my $str = sprintf("%*.*b", $size, $size, $i);
            my %combination;
            for (my $j = 0 ; $j < $size ; $j++)
            {
                if (substr($str, $j, 1)) { $combination{$list[$j]} = $list[$j]; }
            }
            if (scalar keys %combination > 0 && scalar keys %combination <= $size)
            {
                my $flag = 0;
                foreach my $item (@bad)
                {
                    my $flag2 = scalar @{$item};
                    foreach my $key (@{$item})
                    {
                        if (exists $combination{$key}) { $flag2--; }
                    }
                    $flag++ if (!$flag2);
                }
                next if ($flag);
                if ($verbose)
                {
                    foreach my $item (@list) { print {*STDERR} ' ' . $item . ':' . (exists $combination{$item} ? '1' : '0'); }
                    print {*STDERR} "\n";
                }
                my ($fd, $tmp_fasta) = tempfile(UNLINK => 1, OPEN => 0, SUFFIX => '.fasta');
                my $seq_in  = Bio::SeqIO->new(-file => '<' . $fasta,     -format => 'fasta');
                my $seq_out = Bio::SeqIO->new(-file => '>' . $tmp_fasta, -format => 'fasta');
                while (my $inseq = $seq_in->next_seq)
                {
                    my $id = $inseq->id();
                    $id = $1 if ($id =~ m/(^\w+\d+\|\d{4})/);
                    if (exists $sample{$id} && exists $combination{$sample{$id}}) { $seq_out->write_seq($inseq); }
                }
                $seq_out->close();
                $seq_in->close();
                my ($fd2, $tmp_output) = tempfile(UNLINK => 1, OPEN => 0);

                print {*STDERR} $pathlava, '/lava.pl --alignment_fasta ', $tmp_fasta, ' --output_file ', $tmp_output, ' --option_file ', $pathlava, '/t_data/' . ($loose?'loose':'normal') . '_parameters.xml', "\n" if($verbose);
                system($pathlava . '/lava.pl --alignment_fasta ' . $tmp_fasta . ' --output_file ' . $tmp_output . ' --option_file ' . $pathlava . '/t_data/' . ($loose ? 'loose' : 'normal') . '_parameters.xml >/dev/null 2>/dev/null');
                print {*STDOUT} '> ', join(' ', sort keys %combination), "\t";
                if (-r $tmp_output)
                {
                    if (open(my $FILE, $tmp_output))
                    {
                        my $lines = 0;
                        $lines++ while (<$FILE>);
                        close $FILE;
                        print {*STDOUT} int($lines / 12);
                    }
                    print {*STDOUT} "\tPrimer(s) generated!\n";
                    unlink $tmp_output;
                    system('cp ' . $tmp_output . '.dash ' . $fasta . q{_} . join(q{_}, sort keys %combination) . '.primers');
                    unlink $tmp_output . '.dash';
                }
                else
                {
                    my @tmp;
                    foreach my $item (@list) { push @tmp, $item if (exists $combination{$item}); }
                    push @bad, [@tmp];
                    print {*STDOUT} "0\tNo primer generated! Will skip " . join(' ', @tmp) . " in the next steps\n";
                }
                unlink $tmp_fasta;
            }
        }
    }
}
