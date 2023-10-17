sed "s/neutron/gamma/g" LXe.mac >vary_ene/LXe_gamma_1000.mac
./LXe vary_ene/LXe_gamma_1000.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_gamma_1000.root
sed "s/1000 keV/100 keV/g" LXe.mac >vary_ene/LXe_100.mac
./LXe vary_ene/LXe_100.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_100.root

sed "s/1000 keV/500 keV/g" LXe.mac >vary_ene/LXe_500.mac
./LXe vary_ene/LXe_500.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_500.root

sed "s/1000 keV/1000 keV/g" LXe.mac >vary_ene/LXe_1000.mac
./LXe vary_ene/LXe_1000.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_1000.root

sed "s/1000 keV/1500 keV/g" LXe.mac >vary_ene/LXe_1500.mac
./LXe vary_ene/LXe_1500.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_1500.root

sed "s/1000 keV/2000 keV/g" LXe.mac >vary_ene/LXe_2000.mac
./LXe vary_ene/LXe_2000.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_2000.root

sed "s/1000 keV/2000 keV/g" LXe.mac >vary_ene/LXe_2000.mac
./LXe vary_ene/LXe_2000.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_2000.root

sed "s/1000 keV/2500 keV/g" LXe.mac >vary_ene/LXe_2500.mac
./LXe vary_ene/LXe_2500.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_2500.root

sed "s/1000 keV/3000 keV/g" LXe.mac >vary_ene/LXe_3000.mac
./LXe vary_ene/LXe_3000.mac
./addhist.sh
cp LXe_all.root vary_ene/LXe_3000.root
