#!/usr/bin/perl -w

# Copyright (c) 2013
# by Christian Panse <cp@fgcz.ethz.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 NAME 

protViz_mgf2RData.pl - mascot generic  file to RData exporter

=head1 SYNOPSIS

mgf2RData.pl [-m=<mascot generic file>] 

=head1 COPYRIGHT

Copyright (c) 2013 Christian Panse.

GNU General Public License Version 3

=head1 AUTHORS

Christian Panse <cp@fgcz.ethz.ch>

=head1 DESCRIPTION

The program exports Matrix Science (http://www.matrixscience.com/) 
mascot dat files to an R (http://www.r-project.org/) object.

The program reqires R install and is testet on a debian linux system.

The program is part of the protViz R package on CRAN

=head1 OPTIONS

=head2 -m

xxx

=cut

use strict;
use warnings;

sub main(){
    
    my $Rdata = "mgf";
    $Rdata = shift || die "no data name provided";
    my $i = 1;
    my ($_title, $_pepmass, $_charge, $_scan, $_rtinseconds, @mZ, @intensity);


    open (RFILE, " | tee /tmp/dump.R | R --no-save") || die "could not open file for writing ...$!\n";

    print RFILE $Rdata." <- list()\n";

    while (<>){
        s/\r\n/\n/;
        s/\\/\//g;
        chomp;

        if (/^BEGIN IONS/){
            $_scan='NA';
            $_title=$i;
            $_rtinseconds='NA';
        }elsif (/^END IONS/){
            $_charge = "NA";

            print RFILE "\n". $Rdata . "[[" . $i . "]] <- list(\n";
            print RFILE "title=\"" . $_title ."\",\n";
            print RFILE "rtinseconds=" . $_rtinseconds.",\n";
            print RFILE "charge=" . $_charge.",\n";
            print RFILE "scan=" . $_scan.",\n";
            print RFILE "pepmass=" . $_scan.",\n";

            print RFILE "mZ=c(";
            for (my $ii=0; $ii < $#mZ; $ii++){
                print RFILE $mZ[$ii];                    
                print RFILE ", " if ($ii < $#mZ-1);
            }
            print RFILE "),\n";

            print RFILE "intensity=c(";
            for (my $ii=0; $ii < $#intensity; $ii++){
                print RFILE $intensity[$ii];                    
                print RFILE ", " if ($ii < $#intensity-1);
            }
            print RFILE ")\n";

            print RFILE ")\n";
            undef $_title;
            undef $_scan;
            undef $_charge;
            undef $_rtinseconds;
            undef $_pepmass;
            undef @mZ;
            undef @intensity;
            $i = $i +1;

        }elsif (/^TITLE=(.+)/){
            $_title = $1;
        }elsif (/^PEPMASS=(.+)\ (.+)/){
            $_pepmass = $1;
        }elsif (/^CHARGE=(\d+)./){
            $_charge = $1;
        }elsif (/^SCANS=(\d+)/){
            $_scan = $1;
        }elsif (/^RTINSECONDS=(.+)/){
            $_rtinseconds = $1;
        }elsif (/^(\d+\.\d+)\s(\d+\.{0,1}\d+)/){
            push @mZ, $1;
            push @intensity, $2;
        }
   }

   print RFILE "save($Rdata , file='" . $Rdata . ".RData', compress=TRUE)\n";
   close(RFILE);
}

## MAIN

my $a;
while ($a = shift @ARGV) {
    chomp $a;
    if ($a =~ /-n=(.+)/ || $a =~ /--name=(.+)/) {
        &main($1);
    }
}
