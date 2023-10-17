for i in $(seq 200 200 2000)
do
    sed "s/1000 keV/$i keV/g" gmono.mac >g_vary_e/LXe_gamma_$i.mac
    ./LXe g_vary_e/LXe_gamma_$i.mac
    ./haddall.sh LXe_tall.root 24
    cp -rf LXe_tall.root g_vary_e/LXe_gamma_$i.root
done

# sed 's/1000 keV/200 keV/g' gmono.mac >g_vary_e/LXe_gamma_200.mac
# ./LXe g_vary_e/LXe_gamma_200.mac
# ./haddall.sh LXe_tall.root 24
# cp LXe_tall.root g_vary_e/LXe_gamma_200.root

