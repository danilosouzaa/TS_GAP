use strict;
use warnings;
use 5.010;

my @names =(" c05100 100000 5   1931 " , " d05200 100000 5   12742 " , " e10100 100000 5   11577 " , " c10200 100000 5   2806 " , " d10400 100000 5   24963 " , " e15900 100000 5   102421 " , " c20100 100000 5   1243 " , " d20200 100000 5   12241 " , " e20400 100000 5   44877 " , " c30900 100000 5   9982 " , " d40400 100000 5   24392 " , " e60900 100000 5   100153 " , " c201600 100000 5   18802 " , " c401600 100000 5   17145 " , " c801600 100000 5   16285 ");
my $n_lc = 0;
my $n_ex = 0;
foreach my $n (@names){
                open(FIN,">>Res150");
                print FIN ("exp. $n \n");
                close(FIN);
                system("./m2 $n  >>Res150");
}

exit;
