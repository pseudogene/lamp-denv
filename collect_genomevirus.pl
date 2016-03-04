#!/usr/bin/perl
# $Revision: 0.3 $
# $Date: 2016/03/04 $
# $Id: collect_genomevirus.pl $
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
use Bio::SeqIO;
use Bio::DB::GenBank;
my %query = (1 => 'txid11053[Organism:exp] AND "Dengue virus 1"[porgn]', 2 => 'txid11060[Organism:exp] AND "Dengue virus 2"[porgn]', 3 => 'txid11069[Organism:exp] AND "Dengue virus 3"[porgn]', 4 => 'txid11070[Organism:exp] AND "Dengue virus 4"[porgn]');
my ($verbose, $after, $before, $ok, $skip, $missing, $denv, $genome) = (0, 0, 1900 + (localtime(time))[5], 0, 0, 0);
GetOptions('f|fasta=s' => \$genome, 'a|after|y|year:i' => \$after, 'b|before:i' => \$before, 'd|denv=i' => \$denv, 'v|verbose!' => \$verbose);

if (defined $genome && !-e $genome && defined $after && defined $before && defined $denv && int($denv) > 0 && int($denv) < 5 && open(my $seq_out, '>', $genome))
{
    my $gb       = new Bio::DB::GenBank;
    my $thequery = Bio::DB::Query::GenBank->new(-query => $query{int($denv)} . ' AND (complete[All Fields] AND genome[All Fields])' . ($after > 0 ? ' AND ("' . $after . '"[PDAT] : "' . $before . '"[PDAT])' : q{}), -db => 'nucleotide');
    my $seqio    = $gb->get_Stream_by_query($thequery);
    while (my $seq = $seqio->next_seq)
    {
        print {*STDERR} '>', $seq->id(), "\n" if ($verbose);
        my ($date, $country);
        for my $feat_object ($seq->get_SeqFeatures)
        {
            if ($feat_object->primary_tag eq 'source')
            {
                for my $tag ($feat_object->get_all_tags)
                {
                    if ($tag eq 'country')
                    {
                        for my $value ($feat_object->get_tag_values($tag)) { $country = $value; }
                    }
                    elsif ($tag eq 'collection_date')
                    {
                        for my $value ($feat_object->get_tag_values($tag)) { $date = $1 if ($value =~ m/(\d{4})/); }
                    }
                }
            }
        }
        if (defined $date && defined $country)
        {
            if ($date > $after && $date <= $before)
            {
                $ok++;
                print {$seq_out} q{>} . $seq->accession_number() . q{|} . $date . q{|} . $country . "\n" . $seq->seq() . "\n";
            }
            else { $skip++; }
        }
        else { print {*STDERR} '  !No metadata for accession: ', $seq->accession_number(), "\n"; $missing++; }
    }
    close $seq_out;
    print {*STDERR} "\nNumber of sequence retrieved         ", ($ok + $skip + $missing), ($after > 0 ? "\nNumber of sequence too old too young " . $skip : q{}), "\nNumber of sequence without data      ", $missing, "\nNumber of remaining sequences        ", $ok, "\n";
    if (-r $genome)
    {
        system('GramAlign -q -F 1 -f 2 -i ' . $genome . ' -o ' . $genome . '.align.fa');
        unlink('_ga_temp.page0');
        if ($verbose)
        {
            my $seqio = Bio::SeqIO->new('-format' => 'fasta', -file => $genome . '.align.fa');
            my $seqobj = $seqio->next_seq();
            print {*STDERR} 'Alignment length                 ', $seqobj->length, "\n";
        }
    }
}
