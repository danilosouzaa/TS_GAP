use strict;
use warnings;
use 5.010;

my @names =( " e15900 100000 20   102421 " , " c20100 100000 20   1243 " , " d20200 100000 20   12241 " , " e20400 100000 20   44877 " ,  " c201600 100000 20   18802 ");
my $n_lc = 0;
my $n_ex = 0;
foreach my $n (@names){
                open(FIN,">>Resultado20-2");
                print FIN ("exp. $n \n");
                close(FIN);
                system("./new $n  >>Resultado20-2");
}

exit;

