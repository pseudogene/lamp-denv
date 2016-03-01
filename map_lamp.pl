#!/usr/bin/perl
# $Revision: 0.4 $
# $Date: 2016/02/29 $
# $Id: map_lamp.pl $
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
use Bio::Graphics;
use Bio::AlignIO;
use Bio::SeqFeature::Generic;

#----------------------------------------------------------
our ($VERSION) = 0.4;

#----------------------------------------------------------
my ($verbose, $png, $svg, $align) = (0, 0, 0);
my @primerfiles;
GetOptions('p|primer=s' => \@primerfiles, 'a|align|alignment=s' => \$align, 'svg!' => \$svg, 'png!' => \$png, 'v|verbose!' => \$verbose);
if (defined $align && -r $align && scalar @primerfiles > 0)
{
    my $inalign = Bio::AlignIO->new(-file => $align, -format => 'fasta');
    my $aln = $inalign->next_aln() or die 'No alignment';
    my $panel;
    if ($svg) { $panel = Bio::Graphics::Panel->new(-length => $aln->length() * 1.1, -key_style => 'between', -width => 800, -pad_left => 10, -pad_right => 10, -image_class => 'GD::SVG',); }
    else      { $panel = Bio::Graphics::Panel->new(-length => $aln->length() * 1.1, -key_style => 'between', -width => 800, -pad_left => 10, -pad_right => 10,); }
    my $full_length = Bio::SeqFeature::Generic->new(-start => 1, -end => $aln->length());
    $panel->add_track($full_length, -glyph => 'arrow', -bump => 0, -tick => 2, -double => 1, -label => 1,);
    print {*STDERR} '>Alignment length: ', $aln->length(), "\n" if ($verbose);
    my @colors = qw(cyan orange blue purple green chartreuse magenta yellow aqua);
    my $idx    = 0;

    foreach my $pfile (@primerfiles)
    {
        if (open my $IN, '<', $pfile)
        {
            my $track = $panel->add_track(
                -glyph       => 'generic',
                -label       => 1,
                -bgcolor     => $colors[$idx++ % @colors],
                -min_score   => 0,
                -max_score   => 1000,
                -font2color  => 'red',
                -key         => "$pfile",
                -description => sub {
                    my $feature = shift;
                    my $score   = $feature->score;
                    return "penality=$score";
                },
            );
            while (<$IN>)
            {
                if (m/^>(\d+) F3 \(penalty: (\d+).*locations: \((\d+)\-\d+, \d+\-(\d+)\)/)
                {
                    print {*STDERR} $1, q{ }, $3, q{-}, $4, "\n";
                    my $feature = Bio::SeqFeature::Generic->new(-score => $2, -display_name => 'Set ' . ($1 + 1), -start => $3, -end => $4,);
                    $track->add_feature($feature);
                }
            }
        }
    }
    if    ($svg) { print $panel->svg; }
    elsif ($png) { print $panel->png; }
}
