use strict;
use warnings;
use 5.010;

my @names =(" c05100 100000 20   1931 " , " d05200 100000 20   12742 " , " e10100 100000 20   11577 " , " c10200 100000 20   2806 " , " d10400 100000 20   24963 " );
my $n_lc = 0;
my $n_ex = 0;
foreach my $n (@names){
                open(FIN,">>Resultado20-1");
                print FIN ("exp. $n \n");
                close(FIN);
                system("./new $n  >>Resultado20-1");
}

exit;

