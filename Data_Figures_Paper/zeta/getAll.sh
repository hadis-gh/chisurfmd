#! /bin/baish


bash getZeta.sh                ../E_all_SS.dat | tee E_zeta_Rbest_SS.dat
bash getZetaBest.sh        E_zeta_Rbest_SS.dat | tee E_zetaBest_Rbest_SS.dat


bash getZeta.sh                ../E_all_DD.dat | tee E_zeta_Rbest_DD.dat
bash getZetaBest.sh        E_zeta_Rbest_DD.dat | tee E_zetaBest_Rbest_DD.dat
